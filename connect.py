import sys, paramiko

hostname = "10.53.139.90"
password = 'enkf2018'
username = 'dev'

client = paramiko.SSHClient()
client.load_system_host_keys()
client.connect(hostname, username=username, password=password)

stdin, stdout, strerr = client.exec_command("ls ")

datapath = "./wrf/v7.3.1.kc/wrf/20190709/" + which_data

files = 

client.close()


