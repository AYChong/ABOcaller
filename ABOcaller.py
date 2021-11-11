# Convert phased haplotypes to ABO blood groups based on variant genotypes at 
# rs8176719 and rs8176747

import argparse
import os
import sys
import pandas
import numpy

# Create reference table for ABO blood types
blood_types = [["OO", "OA", "AO", "AA"],
				["OO", "OB", "AO", "AB"], 
				["OO", "OA", "BO", "BA"],
				["OO", "OB", "BO", "BB"]]
o_variant = ["T:T", "T:TC", "TC:T", "TC:TC"]
ab_variant = ["C:C", "C:G", "G:C", "G:G"]
abo_df = pandas.DataFrame(blood_types, index=ab_variant, columns=o_variant)

# Dictionaries to convert ABO haplotype/Secretor genotype to simple status
abo_dict = {
	"AA":"A", 
	"AB":"AB", 
	"BA":"AB", 
	"BB":"B", 
	"OA":"A", 
	"OB":"B", 
	"OO":"O", 
	"AO":"A", 
	"BO":"B"
}
se_dict = {"GG":0, "GA":0, "AG":0, "AA":1, "NA":"NA"}

def linecount(sample_file):
	f = open(sample_file, "rb")
	lines = 0
	buf_size = 1024 * 1024
	read_f = f.raw.read
	buf = read_f(buf_size)
	while buf:
		lines += buf.count(b"\n")
		buf = read_f(buf_size)
	return lines

def check_input(infile, **kwargs):
	first_line = ""
	with open(infile, "r") as datafile:
		first_line = datafile.readline()
	if first_line.startswith("##fileformat=VCF"):
		filetype = "vcf"
	else: # File should be oxford format gen/haps file
		sample_file = kwargs.get("sample", None)
		try: 
			# adjust for header and info line in oxford format sample file
			num_samples = linecount(sample_file) - 2 
			gen_length = 5 + (3*num_samples)
			hap_length = 5 + (2*num_samples)
			num_cols = len(first_line.strip().split())
			if num_cols == gen_length:
				filetype = "gen"
			elif num_cols == hap_length:
				filetype = "haps"
			else:
				filetype = "unknown"
		except TypeError:
			print("ERROR:", "check_input:", "Sample file must be specified")
			filetype = "unknown"
	return filetype # one of "gen", "haps", "vcf"

def get_samples_gen(sample_file):
	# read in gen and sample files for variant
	headers = ["chrom", "rsid", "pos", "gen1", "gen2"]
	sample_list = []
	with open(sample_file, "r") as samples:
		for i in range(2):
			next(samples)
		for line in samples:
			fid, iid, _ = line.strip().split(None, 2)
			sid = fid + "_" + iid
			hap1 = sid + "_1"
			hap2 = sid + "_2"
			hap3 = sid + "_3"
			headers.extend([hap1, hap2, hap3])
			sample_list.append(sid)
	return headers, sample_list

def get_samples_vcf(vcf_file):
	sample_list = []
	with open(vcf_file, "r") as samples:
		for line in samples:
			if line.startswith("#CHROM"):
				headers = line.strip().split("\t")
				sample_list = headers[9:]
				break
	return headers, sample_list

def get_abo_type(row, abo_df):
	gen1 = row["O_genotype"]
	gen2 = row["AB_genotype"]
	try:
		abo = abo_df.loc[gen2, gen1]
	except KeyError:
		abo = "NA"
	except Exception as e:
		print(row)
		print(e)
		sys.exit(1)
	return abo

def parse_abo_haps(haps_file, sample_file, o_variant, ab_variant):
	sample_list = []
	o_variant_0 = ""
	o_variant_1 = ""
	ab_variant_0 = ""
	ab_variant_1 = ""
	variant_checklist = []
	
	with open(sample_file, "r") as samples:
		for i in range(2): ## skip double header row in sample file
			next(samples)
		for line in samples:
			fid, iid, _ = line.strip().split(None, 2)
			sid = fid + "_" + iid
			sample_list.append(sid)
	with open(haps_file, "r") as infile:
		for line in infile:
			chrom, variantid, pos, gen1, gen2, indivs = line.strip().split(" ", 5)
			if variantid == o_variant:
				variant_checklist.append(variantid)
				o_variant_0 = gen1
				o_variant_1 = gen2
				hap_line = [numpy.nan if x == "" else int(x) for x in indivs.split(" ")]
				o_variant_0_data = hap_line[::2]
				o_variant_1_data = hap_line[1::2]
			elif variantid == ab_variant:
				variant_checklist.append(variantid)
				ab_variant_0 = gen1
				ab_variant_1 = gen2
				hap_line = [numpy.nan if x == "" else int(x) for x in indivs.split(" ")]
				ab_variant_0_data = hap_line[::2]
				ab_variant_1_data = hap_line[1::2]
			else:
				pass
	if len(variant_checklist) < 1:
		print("ERROR: {0} and {1} not present in haps file!".format(o_variant, ab_variant))
		sys.exit(1)
	else:
		for variant in [o_variant, ab_variant]:
			if not variant in variant_checklist:
				print("ERROR: {0} not present in haps file!".format(variant))
				sys.exit(1)
				
	haplotype_df = pandas.DataFrame({
		"SID":sample_list, 
		"Ovariant_0":o_variant_0_data, 
		"Ovariant_1":o_variant_1_data, 
		"ABvariant_0":ab_variant_0_data, 
		"ABvariant_1":ab_variant_1_data})

	return sample_list, haplotype_df, o_variant_0, o_variant_1, ab_variant_0, ab_variant_1

def split_GT(indivs, gt_idx, gp_idx, gp_threshold, phased):
	hap1_data = []
	hap2_data = []
	indiv_geno_data = indivs.split()
	for geno_data in indiv_geno_data:
		if ":" in geno_data:
			gt_str = geno_data.split(":")[gt_idx]
		else:
			gt_str = geno_data
		if gp_idx:
			gp_str = geno_data.split(":")[gp_idx]
			max_gp = max([float(x) for x in gp_str.split(",")])
		else:
			max_gp = 1 # Assume already thresholded or directly typed
		if max_gp >= gp_threshold and not "." in gt_str:
			if phased:
				try:
					hapA, hapB = gt_str.split("|")
				except ValueError:
					print("ERROR: Data does not appear to be phased... ", gt_str)
			else:
				hapA, hapB = gt_str.split("/")
			hap1_data.append(int(hapA))
			hap2_data.append(int(hapB))
		else:
			hap1_data.append(numpy.nan)
			hap2_data.append(numpy.nan)
	
	return hap1_data, hap2_data

def parse_abo_vcf(vcf_file, o_variant, ab_variant, gp_threshold, phased):
	headers, sample_list = get_samples_vcf(vcf_file)
	
	o_variant_0 = ""
	o_variant_1 = ""
	ab_variant_0 = ""
	ab_variant_1 = ""
	variant_checklist = []
	with open(vcf_file, "r") as infile:
		for line in infile:
			if not line.startswith("#"): #skip info
				chrom, pos, variantid, ref, alt, qual, filt, variantinfo, variantformat, indivs = line.strip().split("\t", 9)
				if ":" in variantformat:
					gt_idx = variantformat.split(":").index("GT")
					if "GP" in variantformat:
						gp_idx = variantformat.split(":").index("GP")
				else: # There should only be GT in the genotype data
					gt_idx = 0
					gp_idx = False
				
				if variantid == o_variant:
					variant_checklist.append(variantid)
					o_variant_0 = ref
					o_variant_1 = alt
					o_variant_0_data, o_variant_1_data = split_GT(indivs, gt_idx, gp_idx, gp_threshold, phased)
				elif variantid == ab_variant:
					variant_checklist.append(variantid)
					ab_variant_0 = ref
					ab_variant_1 = alt
					ab_variant_0_data, ab_variant_1_data = split_GT(indivs, gt_idx, gp_idx, gp_threshold, phased)
				else:
					pass
	if len(variant_checklist) < 1:
		print("ERROR: {0} and {1} not present in vcf file!".format(o_variant, ab_variant))
		sys.exit(1)
	else:
		for variant in [o_variant, ab_variant]:
			if not variant in variant_checklist:
				print("ERROR: {0} not present in vcf file!".format(variant))
				sys.exit(1)
	haplotype_df = pandas.DataFrame({
		"SID":sample_list, 
		"Ovariant_0":o_variant_0_data, 
		"Ovariant_1":o_variant_1_data, 
		"ABvariant_0":ab_variant_0_data, 
		"ABvariant_1":ab_variant_1_data})

	return sample_list, haplotype_df, o_variant_0, o_variant_1, ab_variant_0, ab_variant_1

def call_ABO(infile, outfile, o_variant, ab_variant, gp_threshold, phased, **kwargs):
	sample_file = kwargs.get("sample", None)
	filetype = check_input(infile, sample=sample_file)
	print("call_ABO:", "filetype:", filetype)
	if filetype == "haps":
		sample_list, haplotype_df, o_variant_0, o_variant_1, ab_variant_0, ab_variant_1 = parse_abo_haps(infile, sample_file, o_variant, ab_variant)
	elif filetype == "vcf":
		sample_list, haplotype_df, o_variant_0, o_variant_1, ab_variant_0, ab_variant_1 = parse_abo_vcf(infile, o_variant, ab_variant, gp_threshold, phased)
	else:
		print("ERROR: Input filetype not recognised")
		sys.exit()

	type_df = pandas.DataFrame()
	type_df["SID"] = sample_list
	type_df["Ovariant_0"] = numpy.where(haplotype_df["Ovariant_0"].isnull(), ".", 
										numpy.where(haplotype_df["Ovariant_0"] == 0, o_variant_0, o_variant_1))
	type_df["Ovariant_1"] = numpy.where(haplotype_df["Ovariant_1"].isnull(), ".", 
										numpy.where(haplotype_df["Ovariant_1"] == 0, o_variant_0, o_variant_1))
	type_df["ABvariant_0"] = numpy.where(haplotype_df["ABvariant_0"].isnull(), ".", 
										numpy.where(haplotype_df["ABvariant_0"] == 0, ab_variant_0, ab_variant_1))
	type_df["ABvariant_1"] = numpy.where(haplotype_df["ABvariant_1"].isnull(), ".", 
										numpy.where(haplotype_df["ABvariant_1"] == 0, ab_variant_0, ab_variant_1))
	type_df["O_genotype"] = type_df[["Ovariant_0", "Ovariant_1"]].apply(":".join, axis=1)
	type_df["AB_genotype"] = type_df[["ABvariant_0", "ABvariant_1"]].apply(":".join, axis=1)
	type_df["ABO_haplotype"] = type_df.apply(get_abo_type, args=(abo_df,), axis=1)
	type_df["ABO_bloodtype"] = type_df["ABO_haplotype"].map(abo_dict)
	type_df[["SID", "ABO_haplotype", "ABO_bloodtype", "O_genotype", "AB_genotype"]].to_csv(outfile, sep="\t", header=True, index=False)
	
def parse_se_gen(infile, sample_file, se_variant, gp_threshold):
	headers, sample_list = get_samples_gen(sample_file)
	gen_df = pandas.read_csv(infile, delim_whitespace=True, header=None, names=headers)
	# Check that the correct variant is present
	variant_df = gen_df.loc[gen_df["rsid"] == se_variant]
	if variant_df.empty:
		print("ERROR: {0} not present in gen file!".format(se_variant))
		sys.exit(1)
	
	se_variant_0 = variant_df.loc[variant_df["rsid"] == se_variant, "gen1"].values[0]
	se_variant_1 = variant_df.loc[variant_df["rsid"] == se_variant, "gen2"].values[0]
	
	genotype_df = variant_df.set_index("rsid").drop(["chrom","pos", "gen1", "gen2"], axis=1)
	geno_list = genotype_df.values.tolist()[0]
	geno_tuples = [geno_list[i:i + 3] for i in range(0, len(geno_list), 3)]
	
	genotypes = []
	counter = 0
	for ind in geno_tuples:
		sid = sample_list[counter]
		if max(ind) >= gp_threshold: # Imputed genotype is of sufficient quality
			geno_index = ind.index(max(ind))
			if geno_index == 0: # Homozygous for ref allele
				genotype = se_variant_0 + se_variant_0
			elif geno_index == 1: # Heterozygote
				genotype = se_variant_0 + se_variant_1
			elif geno_index == 2: # Homozygous for alt allele
				genotype = se_variant_1 + se_variant_1
			else:
				genotype = "NA"
		else:
			genotype = "NA"
		genotypes.append(genotype)
		counter += 1
	return sample_list, genotypes, se_variant_0, se_variant_1

def parse_se_vcf(vcf_file, se_variant, gp_threshold, phased):
	headers, sample_list = get_samples_vcf(vcf_file)
	se_variant_0 = ""
	se_variant_1 = ""

	variant_checklist = []
	with open(vcf_file, "r") as infile:
		for line in infile:
			if not line.startswith("#"): # skip info
				chrom, pos, variantid, ref, alt, qual, filt, variantinfo, variantformat, indivs = line.strip().split("\t", 9)
				if ":" in variantformat:
					gt_idx = variantformat.split(":").index("GT")
					if "GP" in variantformat:
						gp_idx = variantformat.split(":").index("GP")
				else: # There should only be GT in the genotype data
					gt_idx = 0
					gp_idx = False
				if variantid == se_variant:
					variant_checklist.append(variantid)
					se_variant_0 = ref
					se_variant_1 = alt
					try:
						se_variant_0_data, se_variant_1_data = split_GT(indivs, gt_idx, gp_idx, gp_threshold, phased)
					except:
						print("ERROR: Are you sure that the data is in vcf format? First genotype is:", indivs.split("\t", 1)[0])
						sys.exit(1)
				else:
					pass
	if len(variant_checklist) < 1:
		print("ERROR: {0} not present in vcf file!".format(se_variant))
		sys.exit(1)
	geno_tuples = zip(se_variant_0_data, se_variant_1_data)
	
	genotypes = []
	counter = 0
	for genA, genB in geno_tuples:
		sid = sample_list[counter]
		try:
			geno_index = int(genA) + int(genB)
			if geno_index == 0: # Homozygous for ref allele
				genotype = se_variant_0 + se_variant_0
			elif geno_index == 1: # Heterozygote
				genotype = se_variant_0 + se_variant_1
			elif geno_index == 2: # Homozygous for alt allele
				genotype = se_variant_1 + se_variant_1
		except ValueError:
			genotype = "NA"
		except Exception as e:
			print("ERROR: Unknown genotype encoding. Treating as null value", sid, genA, genB)
			genotype = "NA"
			print(e)
		genotypes.append(genotype)
		counter += 1

	return sample_list, genotypes, se_variant_0, se_variant_1

def call_Se(infile, outfile, se_variant, gp_threshold, phased, **kwargs):
	sample_file = kwargs.get("sample", None)
	filetype = check_input(infile, sample=sample_file)
	if filetype == "gen":
		sample_list, genotypes, se_variant_0, se_variant_1 = parse_se_gen(infile, sample_file, se_variant, gp_threshold)
	elif filetype == "vcf":
		sample_list, genotypes, se_variant_0, se_variant_1 = parse_se_vcf(infile, se_variant, gp_threshold, phased)
	else:
		print("ERROR: Input filetype not recognised")
	
	with open(outfile, "w") as out:
		out.write("SID" + "\t" + "Se_genotype" + "\t" + "Se_status" + "\n")
		counter = 0
		for sid, genotype in zip(sample_list, genotypes):
			se_factor = se_dict[genotype]
			out.write(sid + "\t" + genotype + "\t" + str(se_factor) + "\n")
			counter += 1


def main():
	#### set up options with parser
	parser = argparse.ArgumentParser()
	parser._action_groups.pop()
	required = parser.add_argument_group("required arguments")
	optional = parser.add_argument_group("optional arguments")
	
	required.add_argument("-i", "--genotype_file", dest="infile", required=True,
						help="File containing genotypes for either ABO or FUT2 variants")
	optional.add_argument("-s", "--sample_file", dest="sample", required=False,
						help="File containing sample IDs. Required if providing .hap or .gen input")
	required.add_argument("-o", "--output_file", dest="outfile", required=True,
						help="File name for the resulting blood type calls")
	required.add_argument("-m", "--mode", dest="mode", choices=["ABO", "Se"], required=True,
						help="Type of blood group to call. Valid options are ABO and Se")
	optional.add_argument("--rs8176719", dest="o_variant", required=False,
						default="rs8176719",
						help="Alternate variant ID for ABO variant rs8176719. Default is rs8176719.")
	optional.add_argument("--rs8176747", dest="ab_variant", required=False,
						default="rs8176747",
						help="Alternate variant ID for ABO variant rs8176747. Default is rs8176747.")
	optional.add_argument("--rs601338", dest="se_variant", required=False,
						default="rs601338",
						help="Alternate variant ID for FUT2 variant rs601338. Default is rs601338.")
	optional.add_argument("--gp_cutoff", dest="gp_threshold", required=False,
						default=0,
						help="Posterior probability threshold for calling genotypes from imputed data. \nUsed when parsing .vcf files and .gen files. Default is 0.8.")
	optional.add_argument("--unphased", dest="phased", action="store_false",
						help="Use to indicate that input data in vcf format is unphased. Default behaviour is to assume vcf data is phased")
	args = parser.parse_args()
	
	# File list
	infile = args.infile
	if not os.path.isfile(infile):
		print("ERROR: Genotype file {0} does not exist".format(infile))
		sys.exit(1)
	if args.sample:
		sample = args.sample
		if not os.path.isfile(sample):
			print("ERROR: Sample file {0} does not exist".format(sample))
			sys.exit(1)
	else:
		sample = None
	outfile = args.outfile
	
	# configs
	mode = args.mode
	o_variant = args.o_variant
	ab_variant = args.ab_variant
	se_variant = args.se_variant
	gp_threshold = float(args.gp_threshold)
	phased = args.phased
	
	if mode == "ABO":
		call_ABO(infile, outfile, o_variant, ab_variant, gp_threshold, phased, sample=sample)
	elif mode == "Se":
		call_Se(infile, outfile, se_variant, gp_threshold, phased, sample=sample)
	else:
		print("ERROR: Specified mode {0} is not a valid mode. Valid options are: 'ABO' or 'Se'".format(mode))
		sys.exit(1)

###################################################
if __name__ == "__main__":
	main()