#!/usr/local/bin/python3
#! python3
"""
dl_data.py

"""

import os, glob, tarfile, datetime, pdb, fnmatch
import zipfile

__author__ = "Alex Chartier"
__copyright__ = "Copyright 2018, Johns Hopkins University Applied Physics Laboratory"


def dl_data(dl_start, dl_end, dl_dir, server_names, datatype='swarm', dirnames='not_smart'):
    # Download data, uncompress and remove old files
    # dl_days: datetimes that you want data for
    # dl_dir: root directory where the data should go
    dl_days = setDownloadDays(dl_start, dl_end)
    dl = downloader(dirnames=dirnames)
    dl.downloadDays = dl_days
    dl.setOutputDirectory(dl_dir)
    dl.readTextInput(server_names)
    dl.downloadData(datatype=datatype)
    dl.uncompress(dirnames=dirnames)


class downloader(object):
   """ Data downloading object to get all IDA data files"""
   platform = None

   def __init__(self, dirnames, bin_dir='./'):
      self.dirnames = dirnames
      self.setStartDay()
      self.setSysCommands(bin_dir)

   def setSysCommands(self, bin_dir):
      self.SysCommands = {}

      # -r: recursive, -nv: non-verbose, -nd: no directory structure, -N: only download if newer
      # Reject md5 files - don't know what they are, but I don't want them.
      self.SysCommands['wget_arg'] = 'wget --reject md5 -r -nv -nd -N '
      self.SysCommands['crx2rnx'] = bin_dir + '/CRX2RNX '

   def downloadData(self, datatype='GPS'):
      assert hasattr(self,'outputDirectory'), 'outputDirectory not set'
      os.system('mkdir -p ' + self.outputDirectory)
      one_day = datetime.timedelta(days=1)
      days = [self.downloadDays[0] - one_day] + self.downloadDays + [self.downloadDays[-1] + one_day]
      
      if datatype is 'GPS':
          # Get JPLG* and BRDC* first, 
          ephem_key = [s for s in self.servers if 'brdc' in s][0]
          bias_key = [s for s in self.servers if 'bias' in s][0]

          self.downloadServer(self.servers[ephem_key], days)
          self.downloadServer(self.servers[bias_key], days)

          # Copy across the ephemeris files from previous and next days
          for day in self.downloadDays:
             day_before = day - one_day
             day_after = day + one_day
             if self.dirnames == 'smart':
                dir_str = self.outputDirectory + '/%Y-%j/unprocessed/'
             else:
                dir_str = self.outputDirectory
             # Check for existence of jplg and brdc files 
             brdc_str = dir_str + 'brdc%j* '
             jplg_str = dir_str + 'jplg%j* '
             fnames = day_before.strftime(brdc_str), day.strftime(brdc_str), day_after.strftime(brdc_str), \
                      day_before.strftime(jplg_str), day.strftime(jplg_str), day_after.strftime(jplg_str), 
             [checkFileExists(fname) for fname in fnames]

             out1 = os.system('cp ' + day_before.strftime(brdc_str) + day.strftime(dir_str))
             out2 = os.system('cp ' + day_after.strftime(brdc_str) + day.strftime(dir_str))
             assert (out1 == 0) and (out2 == 0), 'copy failed'

      for day in self.downloadDays:
         for name, server in self.servers.items():
            self.downloadServer(server, [day])

   def downloadServer(self, server, days):
      print("Downloading from " + server.Name)
      filesBefore = len(glob.glob(self.outputDirectory.strip() + '*'))
      
      Command = self.SysCommands['wget_arg']
      if server.Username != '':
         Command += ' --user=' + server.Username + ' --password=' + server.Password + ' '

      for date in days:
         for path in server.DataPaths:
            for pattern in server.Pattern:
               Cmd = Command + '-A ' + date.strftime(pattern) + ' -P ' + self.outputDirectory 
               if self.dirnames == 'smart':
                  Cmd += date.strftime('%Y-%j/unprocessed/') 
               Cmd += ' ' + server.Address + date.strftime(path) + '/'
               os.system(Cmd)
               
      server.filesRetrieved = len(glob.glob(self.outputDirectory.strip() + '*')) - filesBefore

   def setStartDay(self, startDay=-14):
      self.startDay = startDay
      Today = datetime.datetime.utcnow()
      Today = Today.replace(hour=0, minute=0, second=0, microsecond=0)
      self.downloadDays = [Today + datetime.timedelta(X) for X in range(self.startDay, 1)]
      return "Download files newer than %s days" % (self.startDay)

   def readTextInput(self, textFile):
      """ Read the list of servers"""
      with open(textFile) as F: Entries = F.readlines()
      Entries = [Entry.strip('\n') for Entry in Entries]  # Remove newlines
      Entries = [Entry.strip(' ') for Entry in Entries]  # Remove whitespace
      Entries = [Entry for Entry in Entries if not Entry.startswith('#')]  # remove comments
      Entries = filter(None, Entries)  # remove empty lines

      for Entry in Entries:
         addServer_Args = Entry.split(',')
         addServer_Args = [e.strip(' ') for e in addServer_Args]
         addServer_Args = [e.replace("'", "") for e in addServer_Args]
         addServer_Args = [e.replace('"', "") for e in addServer_Args]
         downloader.addServer(self, *addServer_Args)

   def addServer(self, Name, Address, DataPaths, Pattern='', Username='', Password=''):
        """ Create data source objects for the servers"""
        DataPaths = DataPaths.split(';')
        Pattern = Pattern.split(';')
        DataPaths = [d.strip() for d in DataPaths]
        Pattern = [d.strip() for d in Pattern]
        Server = data_source(Name, Address, DataPaths, Pattern, Username, Password, self.downloadDays, self.SysCommands)
        if hasattr(self,'servers'): self.servers[Name] = Server
        else: self.servers = {Name : Server}

   def setOutputDirectory(self,outputDirectory='./input_files'):
        self.outputDirectory = outputDirectory
   
   def uncompress(self, dirnames='smart'):
      """Uncompress all the files we have in the directories"""
      """Maintain original modification date"""

      for day in self.downloadDays:
         # Specify download directory
         if dirnames == 'smart':
            dayDir = self.outputDirectory + day.strftime('%Y-%j/unprocessed/')
         else:
            dayDir = self.outputDirectory
         
         dayDir = dayDir.strip()
         print("Decompressing %s" % dayDir)

         # List all files
         zFiles = []
         dFiles = []
         oFiles = []
         gzFiles = []
         zipFiles = []  
         for root, dn, filenames in os.walk(dayDir):
            for filename in fnmatch.filter(filenames, '*.Z'):
               zFiles.append(os.path.join(root, filename))
            for filename in fnmatch.filter(filenames, '*d'):
               dFiles.append(os.path.join(root, filename))
            for filename in fnmatch.filter(filenames, '*o'):
               oFiles.append(os.path.join(root, filename))
            for filename in fnmatch.filter(filenames, '*.gz'):
               gzFiles.append(os.path.join(root, filename))
            for filename in fnmatch.filter(filenames, '*.ZIP'):
               zipFiles.append(os.path.join(root, filename))
            for filename in fnmatch.filter(filenames, '*.zip'):
               zipFiles.append(os.path.join(root, filename))
            # Consider *i and *n files as *d files for this section
            for filename in fnmatch.filter(filenames, '*i'):
               dFiles.append(os.path.join(root, filename))
            for filename in fnmatch.filter(filenames, '*n'):
               dFiles.append(os.path.join(root, filename))

            # .gz files (gzip)
            self.new_gzFiles = []
            for fileName in gzFiles:
               if fileName[0:-3] not in dFiles: self.new_gzFiles.append(fileName)
               print(str(len(self.new_gzFiles)), ".gz files")
         
            for fileName in self.new_gzFiles:
               tFile = tarfile.open(fileName, 'r:gz', errorlevel=1)
               tFile.extractall(dayDir)

            # .zip files
            self.new_zipFiles = []
            for fileName in zipFiles:
               if fileName[0:-4] not in dFiles: self.new_zipFiles.append(fileName)
            print(str(len(self.new_zipFiles)), ".zip files")
            for fileName in self.new_zipFiles:
               zipfile.ZipFile(fileName).extractall(dayDir)
        
            # .Z files (gzip)
            self.new_zFiles = []
            for fileName in zFiles:
               if fileName[0:-2] not in dFiles: self.new_zFiles.append(fileName)
            print(str(len(self.new_zFiles)), "new .Z files")
            for fileName in self.new_zFiles:
               os.system('gzip -fd ' + fileName + ' > /dev/null')

            # d / RINEX files (CRX2RNX)
            dFiles = [] # Now forget about the *i and *n files
            oFiles = []
            print("Finished decompressing. Compiling list of new d files...")
            os.system('find %s -size  0 -print0 |xargs -0 rm' % dayDir)  # Remove empty files first
            for root, dn, filenames in os.walk(dayDir):
               for filename in fnmatch.filter(filenames, '*d'):
                  dFiles.append(os.path.join(root, filename))
               for filename in fnmatch.filter(filenames, '*o'):
                  oFiles.append(os.path.join(root, filename))
            
            self.new_dFiles = list(set([f[0:-1] for f in dFiles]) - set([f[0:-1] for f in oFiles]))
            self.new_dFiles = [f + 'd' for f in self.new_dFiles]
            self.new_dFiles.sort()
            if self.new_dFiles:
               print("Converting", str(len(self.new_dFiles)), "new d files to o files")
               assert os.path.exists(self.SysCommands['crx2rnx'].strip()), \
                  self.SysCommands['crx2rnx'] + 'executable is not present'

            for fileName in self.new_dFiles:
               print(fileName + ' ----> o')
               os.system(self.SysCommands['crx2rnx'] + fileName + '> /dev/null')
               ofileName = os.path.basename(fileName)[:-1] + 'o'
               oDir = os.path.dirname(fileName)
               shutil.move(ofileName, oDir)
               os.popen('touch -r ' + fileName +' '+ fileName[0:-1] + 'o') # Keep timestamps

            self.oFiles = []  
            for root, dn, filenames in os.walk(dayDir):
               for filename in fnmatch.filter(filenames, '*o'):
                  self.oFiles.append(os.path.join(root, filename))
            print("Finished decompression.", len(self.oFiles), "oFiles present")
            print(" ")


   def purge(self, fileTypes, dirnames='smart'):
      if isinstance(fileTypes, str): 
         fileTypes = [fileTypes,]
      assert isinstance(fileTypes, list), "fileTypes must be list"

      # initialise empty lists if necessary
      self.oldFiles = [] if not hasattr(self, 'oldFiles') else None
      self.emptyDirs = [] if not hasattr(self, 'emptyDirs') else None
      
      for day in self.downloadDays:
         # Specify download directory
         if dirnames == 'smart':
            dayDir = self.outputDirectory + day.strftime('%Y-%j/')
         else:
            dayDir = self.outputDirectory

         dayDir = dayDir.strip()
         print("Listing files in %s" %dayDir)
         for fileType in fileTypes:
            for root, dirnames, filenames in os.walk(dayDir):
               for filename in fnmatch.filter(filenames, fileType):
                  self.oldFiles.append(os.path.join(root, filename))

         print('Deleting', str(len(self.oldFiles)))
         [os.remove(file) for file in self.oldFiles]

      Dirs = []
      for dl_dir in os.listdir(self.outputDirectory):
         if os.path.isdir(dl_dir):
            Dirs = os.path.join(self.outputDirectory, dl_dir)
            self.emptyDirs = [dl_dir for dl_dir in Dirs if os.listdir(dl_dir) == []]

      print('Removing ' + str(len(self.emptyDirs)) + ' empty directories')
      [os.removedirs(dl_dir) for dl_dir in self.emptyDirs]

   def report(self):
      """ Report back process statistics """
      print('   ')
      for name, server in self.servers.items():
         if hasattr(server, 'filesRetrieved'):
            print(server.Name + ': ' + str(server.filesRetrieved) + ' files retrieved')
         else:
            print(server.Name + ': no files retrieved')
         print('   ')

         fileTypes = ['new_zipFiles', 'new_gzFiles', 'new_zFiles', 'new_dFiles']
         for type in fileTypes:
            if hasattr(self,type):
               print(str(len(getattr(self,type))) +' '+ type + ' uncompressed')
            else:
               print('No ' + type + ' uncompressed')
         print('   ')

         if hasattr(self,'oldFiles'):
            print(str(len(self.oldFiles)) + ' old files deleted')
         else:
            print('No files deleted')

         if hasattr(self,'emptyDirs'):
            print(str(len(self.emptyDirs)) + ' empty directories removed')
         else:
            print('No directories removed')
         print('   ')


class data_source(object):
   """ Stores server information, e.g. GPS FTP addresses etc."""
   def __init__(self, Name, Address, DataPaths, Pattern, Username, Password, downloadDays, SysCommands):
      self.Name = Name
      self.Address = Address
      self.DataPaths = DataPaths
      self.Pattern = Pattern
      self.Username = Username
      self.Password = Password
      self.startDay = min(downloadDays) - max(downloadDays)
      self.SysCommands = SysCommands

      if type(self.DataPaths) is str: self.DataPaths = [self.DataPaths]
      if type(self.Pattern) is str: self.Pattern = [self.Pattern]


def processGPS(downloadDays, outputDirectory, bin_dir, inp_dir, dirnames='smart'):
  Count = -1
  for date in downloadDays:
     # Define directory structures
     Count += 1
     if dirnames == 'smart':
        dayDir = outputDirectory + date.strftime('%Y-%j/unprocessed/')
        dayDir = dayDir.strip()
        outDir = dayDir + '../processed/'
     else:
        dayDir = outputDirectory
        dayDir = dayDir.strip()
        outDir = dayDir

     if not os.path.exists(outDir):
        os.makedirs(outDir)
   
     # Copy over a JPL file from previous day if necessary
     jplFile = ''.join(glob.glob(dayDir.strip() + 'jplg*i'))
     if not os.path.exists(jplFile):
        previousDay = date - datetime.timedelta(days=1)
        if dirnames == 'smart':
           previous_dayDir = outputDirectory + previousDay.strftime('%Y-%j/unprocessed/') 
        else:
           previous_dayDir = outputDirectory
        
        previous_jplFile = ''.join(glob.glob(previous_dayDir.strip() + 'jplg*i')) 
        if os.path.exists(previous_jplFile):
           os.system(" ".join(['cp', previous_jplFile, dayDir]))
           print('Copied JPL file from previous day')
        else:
           print('No JPL File - cannot process GPS')
     else:
        print('Found a JPL file')
  
     if os.path.exists(jplFile):
         jplProcFile = calculateBias(jplFile)

         # List relevant files
         rinexFiles = glob.glob(dayDir + '*o')
         procFiles = glob.glob(outDir + '*nc')
        
         """
         new_rinexFiles = list(set([os.path.basename(f) for f in rinexFiles]) \
            - set([os.path.splitext(os.path.basename(f))[0] for f in procFiles]))
         new_rinexFiles = [dayDir + f for f in new_rinexFiles]
         new_rinexFiles.sort()
         """
         rinexFiles.sort()
         brdcFile = glob.glob(dayDir + 'brdc*n')
         assert (isinstance(brdcFile, list)) & (not not brdcFile), 'No Bias File - stopping'
         brdcFile = brdcFile[0] 
         
         # Process RINEX files
         for rinexFile in rinexFiles:
            try:
                print(rinexFile)
                with open('proc_errs.txt', 'w') as f:
                    SGL_igs2tecProcess.process(date, dayDir, bin_dir, inp_dir, \
                         rinexFile, useZeroBias=False, trimRINEX=True, fout=f, dirnames=dirnames)
            except:
                print('Failed to convert %s.' % rinexFile)


def calculateBias(jplFile):
  #%% Possibly not needed because GPStk does this?
  # GSB 07/19/17: Yes needed.  And biases in original JPL file are in nanoseconds.
  # we have to convert to TECU.  1 ns = 2.86 TECU
  "Extract bias estimates from JPL files and write out"

  # Read biases from file
  with open(jplFile) as f: jplText = f.readlines()
  
  i = 0
  for line in jplText:
     Entries = line.split()
     if Entries[0]=='DIFFERENTIAL':
        startLine = i + 1
     if len(Entries) > 5 and Entries[3] == 'STATION' and Entries[5] == 'BIAS':
        endLine = i
        break
     i += 1
  Biases = jplText[startLine:endLine]

  # Write biases to new file
  procBiasFile = jplFile + ".proc"
  with open(procBiasFile, 'w') as f:
     for entry in Biases:
        entryStrings = entry.split()
        PRN = entryStrings[0]
        Bias = float(entryStrings[1])*2.86
        f.write("%s %.3f\n" %(PRN, Bias))

  return procBiasFile


def setDownloadDays(startDay, endDay):
  delta = datetime.timedelta(days=1)
  downloadDays = [startDay,] 
  day = startDay
  while day < endDay:
     day += delta
     downloadDays.append(day)
  return downloadDays


def checkFileExists(fname):
    assert len(glob.glob(fname.strip())) > 0, 'File %s does not exist' % fname


if __name__ == '__main__':
    main()
