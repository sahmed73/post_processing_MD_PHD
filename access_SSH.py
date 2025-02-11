# -*- coding: utf-8 -*-
"""
Author: Shihab Ahmed
Created on Thu Sep 12 00:46:49 2024
"""

import paramiko

# Set up SSH client
ssh = paramiko.SSHClient()
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())

# Connect to the HPC via SSH
hostname = 'login.rc.ucmerced.edu'
username = 'sahmed73'
password = '$i1oveMyself$'  # Alternatively, you can use an SSH key
ssh.connect(hostname, username=username, password=password)
#%%

# Use SFTP to transfer files
sftp = ssh.open_sftp()

# Download a file from the HPC to your local machine
remote_file = '/remote/hpc/path/file.csv'
local_file = '/local/path/file.csv'
sftp.get(remote_file, local_file)  # Download file from HPC

# You can also upload files to the HPC
# sftp.put(local_file, remote_file)

# Close the SFTP connection and SSH session
sftp.close()
ssh.close()

print(f"File downloaded from {remote_file} to {local_file}")
