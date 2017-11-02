import os
import argparse	
def vcf_gz_readin(file_in,source_name):
	fin=os.popen(r'''zcat %s '''%(file_in))
	header=[]
	info=[]
	for line in fin:
		pin=line.strip().split()
		if pin[0][:2]=='##': header.append(pin)
		elif pin[0][0]=='#': info.append(pin)
		else:
			if not source_name in pin[2]:
				pin[2]='_'.join([pin[2].split('_')[0],source_name]+pin[2].split('_')[1:])
			info.append(pin)
	fin.close()
	return [header, info]

def vcf_info_write(fileout, header, info):
	fo=open(fileout,'w')
	for i in header: print(' '.join(i), file=fo)
	for i in info:	print('\t'.join(i), file=fo)
	fo.close()

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser = argparse.ArgumentParser("modify_SV_ID.py")
    parser.add_argument("input", type=str, help="namd of input vcf.gz to be modified")
    parser.add_argument("source", type=str, help="name of algorithm that processed the vcf")
    args = parser.parse_args()
    [header, info]=vcf_gz_readin(args.input,args.source)
    vcf_info_write(args.input.replace('.vcf.gz','.SV_ID_Modi.vcf'), header, info)
    os.system(r'''bgzip %s'''%(args.input.replace('.vcf.gz','.SV_ID_Modi.vcf')))
    os.system(r'''tabix %s'''%(args.input.replace('.vcf.gz','.SV_ID_Modi.vcf.gz')))

if __name__ == '__main__':
    main()


