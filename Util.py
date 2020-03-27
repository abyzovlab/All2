import os
import sys
import pwd
import time
import gzip

def FileValidator(file):
	if os.path.exists(file)==False:
		raise Exception("\nERROR:File: "+file+" Doesn't exists")
	if os.path.isfile(file)==False:
                raise Exception( '\nERROR:File expected:'+file+':is not a file')
	if os.stat(file)[6]==0:
		raise Exception( '\nERROR:File:'+file+':is empty')
	else:
		lines=0
		for i in open(file):
			if i.startswith("#"):
				continue
			lines=1
			break
		if lines==0:
			raise Exception( '\nERROR:File:'+file+':is empty')
			
	if os.access(file,os.R_OK)==False:
                raise Exception( '\nERROR:File:\n'+file+':no read access ')
	return file

def DirectoryValidator(file,read_write="NO"):
        if os.path.exists(os.path.dirname(file))==False:
                raise Exception( '\nERROR:Path:'+file+'\n:Does not exist')
        if read_write=="YES":
                if os.access(os.path.dirname(file),os.W_OK)==False:
                        raise Exception('\nERROR:File:\n'+os.path.dirname(file)+':no write access ')
        if os.path.isfile(file)==True:
                raise Exception('\nERROR:Directory expected:\n'+file+':is a file')
        return file

def check_vcf_with_chr_or_not(vcf_file):
        '''checks if a vcf has "chr" before it's chromosome or not'''
        if vcf_file.endswith(".vcf"):
                vcf_file_fh=open(vcf_file)
        elif vcf_file.endswith(".gz"):
                vcf_file_fh=gzip.open(vcf_file)
        for i in vcf_file_fh:
                if i.startswith("#"):
                        continue
                if i.startswith("chr"):
                        vcf="with_chr"
                else:
                        vcf="without_chr"
                break
        vcf_file_fh.close()
        return vcf

def GenerateErrorMsg(script_name,error_msg,std_err):
	msg="Error encountered while running the script : "+script_name+"\n\n"
	envir_file=".".join(std_err.split(".")[:4])+".sge_stats"
	c=open(envir_file,'r')
	for i in c:
		msg=msg+i
	msg=msg+"\\n\n"+str(error_msg)
	msg=msg+"\\n\n"+"Script to check:"+script_name
	msg=msg+"\\n\n"+"Please check the log file and generate/correct the file"
	msg=msg+"\\n\n"+"Once done restart the job by using .qmod -c jobname"
	msg=msg+"\\n\n"+"Please contact Sarangi.vivekananda@mayo.edu for further assistance."
	return msg
	
	
def OsEnviron(envir):
	if "JOB_ID" not in envir:
		return 0
	envir_file=".".join(envir['SGE_STDERR_PATH'].split(".")[:4])+".sge_stats"
	c=open(envir_file,'w')
	for i in ['JOB_ID','JOB_NAME','HOSTNAME','SGE_STDOUT_PATH','SGE_STDERR_PATH','SGE_O_WORKDIR']:
        	c.write(i+"="+envir[i]+"\n")
	c.close()
	
def SendEmail(script_name,error_message,stderr):
	subject="Error in the Low coverage QC pipeline"
	message=GenerateErrorMsg(script_name,error_message,stderr)
	email=pwd.getpwuid(os.getuid())[4] 
	cmd=" ".join(["echo -e", "\""+message+"\"","| /bin/mailx -s","\""+subject+"\"", "\""+email+"\"",";sleep 10"])
	os.popen(cmd)
	time.sleep(10)

	
def ParseConfig(file):
	config={}
	StoreConfig(file,config)
	return config

def ensure_dir(f):
        '''Checks if a directory is present or not and creates one if not present'''
        if not os.path.exists(f):
                os.makedirs(f)

if __name__ == "__main__":
	try:
        	print(ParseConfig(sys.argv[1]))
	except:
		print(sys.exc_info()[1])
		exit()
