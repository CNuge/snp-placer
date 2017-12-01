import unittest

def snp_contig_location(flag, pos, adjusted_bp_location, alignment_length):
	""" determine new bp position of the snp on the larger contig"""
	try:
		adjusted_bp_location / 1
	except:
		return '-' 
	if flag == 0 or flag == 256:
		""" forward aligment, add adj_bp to pos"""
		return (pos + adjusted_bp_location - 1 )
	elif flag == 16 or flag == 272:
		return (pos + alignment_length - adjusted_bp_location)
	else:
		return '-'

def compliment_name(name, flag):
	""" if the alignment is a reverse, add _comp to the end of its identification """
	if flag == 16 or flag == 272:
		return '%s_comp' % (name)
	else:
		return name
		
def match_snp(in_allele):
	if in_allele == 'A':
		return 'T'
	elif in_allele == 'T':
		return 'A'
	elif in_allele == 'C':
		return 'G'
	elif in_allele == 'G':
		return 'C'

def allele_comp_check(in_allele, flag):
	""" if alignment is a reverse, flip the allele to the complimentary base pair. """
	if flag == 0 or flag == 256:
		return in_allele
	elif flag == 16 or flag == 272:
		if len(in_allele) == 1:
			return match_snp(in_allele)
		else:
			in_alleles = in_allele.split(',')
			out_alleles = []
			for i in in_alleles:
				out_alleles.append(match_snp(i))
			return ','.join(out_alleles)




if __name__ == '__main__':
