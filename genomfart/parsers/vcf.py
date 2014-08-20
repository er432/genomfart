import re
import sys
import gzip

if sys.version_info[0] > 2:
    xrange = range

## Define some regular expressions that can be used to recognize certain data types
# Integers
int_re = re.compile(r'^-*[\d]+$')
# Floats
float_re = re.compile(r'(^-*[\d]+\.[\d]+$|^-*[\d]{1}([\.]?[\d])*(E|e)-[\d]+$)')
# Multiple integers
multi_int_re = re.compile(r'^((-*[\d]+),*)+$')
# Multiple floats
multi_float_re = re.compile(r'^((-*[\d]+\.[\d]+$|^-*[\d]{1}([\.]?[\d])*(E|e)-[\d]),*)+$')
## Parser used for VCF (v4.1) files
class VCF_parser:
    """
    Parser for VCF files
    """
    def __init__(self, vcf_file):
        """
        Instantiates a parser for the VCF file

        Parameters
        ----------

        vcf_file : str
            The path to a VCF file (May or may not be gzipped)
        """
        self.version = None
        self.info_dict = {}
        self.header_dict = {}
        self.genotypes = set()
        self.current_line = 0
        self.start_genotype_byte = 0
        if vcf_file.endswith('.gz'):
            self.file_handle = gzip.open(vcf_file)
        else:
            self.file_handle = open(vcf_file)
        # Parse through the beginning of the file
        while 1:
            line = self.file_handle.readline()
            if line.startswith('##fileformat='):
                # Get VCF version
                self.version = line.split('=')[1]
            elif line.startswith('##INFO='):
                line = line.strip()
                # Get field information
                info = '='.join(line.split('=')[1:])[1:-1]
                non_quoted_fields = list(map(lambda x: x.split('='),\
                                        re.findall("(?<=,)[^\s=]+=[^\"\s,=]+",info)))
                quoted_fields = list(map(lambda x: x.split('='),\
                                    re.findall("(?<=,)[^\s=]+=\"[^\"]+\"",info)))
                info_id = re.search("(?<=ID=)[^\s,]+",info).group()
                self.info_dict[info_id] = dict((k,v) for k,v in list(non_quoted_fields+quoted_fields))
            elif line.startswith('#CHROM'):
                # Get the header and the genotypes available
                line = line[1:].strip().split('\t')
                self.header_dict = dict((k,i) for i,k in enumerate(line))
                self.genotypes = line[9:]
                # Set the starting genotype byte
                self.start_genotype_byte = self.file_handle.tell()
                # Break from the loop
                break
    def parse_geno_depths(self):
        """
        Iterates through the VCF file, getting the genotypes at each position

        Returns
        -------
        A generator that generates tuples of chrom,pos,(ref,alt1,alt2,...),
        {sample->(base_depths)}
        """
        # Go to the start of the genotyping part of the file
        file_seek_pos = self.file_handle.seek(self.start_genotype_byte,0)
        # Go through the file
        for line in self.file_handle:
            line = line.strip().split('\t')
            if len(line) < 2:
                continue
            # Get chromosome and position
            chrom,pos = map(int,line[:2])
            # Get the possible alleles
            alt_alleles = line[self.header_dict['ALT']]
            if alt_alleles != '.':
                alleles = tuple([line[self.header_dict['REF']]]+\
                            line[self.header_dict['ALT']].split(','))
            else:
                alleles = (line[self.header_dict['REF']],)
            # Get the index of the allele depth field
            AD_ind = line[self.header_dict['FORMAT']].split(':').index('AD')
            # Get the depths of each genotype
            geno_dict = {}
            for sample in self.genotypes:
                samp_line = line[self.header_dict[sample]]
                if samp_line.startswith('./.'): geno_dict[sample] = tuple(map(lambda x: 0,alleles))
                else:
                    samp_line = samp_line.split(':')
                    geno_dict[sample] = tuple(map(int,samp_line[AD_ind].split(',')[:len(alleles)]))
            yield chrom,pos,alleles,geno_dict
            self.current_line += 1
    def parse_select_geno_depths(self, genos, info_dict = False, use_chrom = None,
                                 start = None, end = None):
        """
        Iterates through the VCF file, getting selected genotypes at each position.
        Note thtat this assumes samples contain depths

        Parameters
        ----------
        genos : list
           The names of the genotypes you want
        info_dict : boolean
            Whether you want the info dict on the end of the return
        use_chrom : str, optional
            Optional chromosome on which to start the scan
        start : int, optional
            Optional place to start the scan. If chrom also specified, it will be
            within the chromsoome. Otherwise, it will be within the first chromosome
            to have at least the start point
        end : int, optional
            Optional place to end the scan. (Nested within chrom if specified)
        
        Returns
        -------
        A generator that generates tuples of chrom,pos,(ref,alt1,alt2,...),{sample->base_depths),
                                            <info_dict if desired>}
        """
        # Go to the start of the genotyping part of the file
        file_seek_pos = self.file_handle.seek(self.start_genotype_byte,0)
        # Go through the file
        for line in self.file_handle:
            line = line.strip().split('\t')
            if len(line) < 2:
                continue
            # Get chromosome and position
            chrom,pos = line[0],int(line[1])
            if (use_chrom is not None and chrom != use_chrom):
                continue
            elif start and not end:
                if pos < start: continue
            elif end and not start:
                if pos > end: continue
            elif start and end:
                if not start <= pos <= end:
                    continue
            # Get the possible alleles
            alt_alleles = line[self.header_dict['ALT']]
            if alt_alleles != '.':
                alleles = tuple([line[self.header_dict['REF']]]+\
                            line[self.header_dict['ALT']].split(','))
            else:
                alleles = (line[self.header_dict['REF']],)
            # Get the index of the allele depth field
            AD_ind = line[self.header_dict['FORMAT']].split(':').index('AD')
            # Get the depths of each genotype
            geno_dict = {}
            for sample in genos:
                samp_line = line[self.header_dict[sample]]
                if samp_line.startswith('./.'): geno_dict[sample] = tuple(map(lambda x: 0,alleles))
                else:
                    samp_line = samp_line.split(':')
                    geno_dict[sample] = tuple(map(int,samp_line[AD_ind].split(',')[:len(alleles)]))
            ## Create dictionary for info if necessary
            line_info_dict = {}
            if info_dict:
                for field in line[self.header_dict['INFO']].split(';'):
                    if '=' in field:
                        key,val = field.split('=')
                        # Check for regular expressions
                        if int_re.match(val):
                            val = int(val)
                        elif float_re.match(val):
                            val = float(val)
                        elif multi_int_re.match(val):
                            val = tuple(map(int, val.split(',')))
                        elif multi_float_re.match(val):
                            val = tuple((map(float,val.split(','))))
                        line_info_dict[key] = val
                    else:
                        key = field
                        line_info_dict[key] = None
            if info_dict:
                yield chrom,pos,alleles,geno_dict,line_info_dict
            else:
                yield chrom,pos,alleles,geno_dict
            self.current_line += 1
    def parse_select_geno_generic(self, genos, info_dict = False, use_chrom = None,
                                   start = None, end = None):
        """
        Iterates through the VCF file, getting selected genotypes at each position.
        This makes no assumptions about the format of the sample information for each genotype

        Parameters
        ----------

        genos : list
            The names of the genotypes you want
        info_dict : dict
            Whether you want the info dict on the end of the return
        use_chrom : str, optional
            Optional chromosome on which to start the scan
        start : int, optional
            Optional place to start the scan
        end : int, optional
            Optional place to end the scan
        
        Returns
        -------
        A generator that generates tuples of chrom,pos,(ref,alt1,alt2,...),{sample->{prefix->val}),
                                            <info_dict if desired>}
        """
        # Go to the start of the genotyping part of the file
        file_seek_pos = self.file_handle.seek(self.start_genotype_byte,0)
        # Go through the file
        for line in self.file_handle:
            line = line.strip().split('\t')
            if len(line) < 2:
                continue
            # Get chromosome and position
            chrom,pos = line[0],int(line[1])
            if (use_chrom is not None and chrom != use_chrom):
                continue
            elif start and not end:
                if pos < start: continue
            elif end and not start:
                if pos > end: continue
            elif start and end:
                if not start <= pos <= end:
                    continue
            # Get the possible alleles
            alt_alleles = line[self.header_dict['ALT']]
            if alt_alleles != '.':
                alleles = tuple([line[self.header_dict['REF']]]+\
                            line[self.header_dict['ALT']].split(','))
            else:
                alleles = (line[self.header_dict['REF']],)
            # Get info of each genotype
            geno_keys = line[self.header_dict['FORMAT']].split(':')
            geno_dict = {}
            for sample in genos:
                samp_line = line[self.header_dict[sample]].split(':')
                try:
                    geno_dict[sample] = dict((k,samp_line[i]) for i,k in enumerate(geno_keys))
                except IndexError:
                    geno_dict[sample] = {geno_keys[0]: samp_line[0]}
            ## Create dictionary for info if necessary
            line_info_dict = {}
            if info_dict:
                for field in line[self.header_dict['INFO']].split(';'):
                    if '=' in field:
                        key,val = field.split('=')
                        # Check for regular expressions
                        if int_re.match(val):
                            val = int(val)
                        elif float_re.match(val):
                            val = float(val)
                        elif multi_int_re.match(val):
                            val = tuple(map(int, val.split(',')))
                        elif multi_float_re.match(val):
                            val = tuple((map(float,val.split(','))))
                        line_info_dict[key] = val
                    else:
                        key = field
                        line_info_dict[key] = None
            if info_dict:
                yield chrom,pos,alleles,geno_dict,line_info_dict
            else:
                yield chrom,pos,alleles,geno_dict
            self.current_line += 1
    def parse_site_infos(self, filter_excludes = None, filter_requires = None):
        """
        Iterates through the VCF file, getting the info for each site

        Parameters
        ----------
        filter_excludes : set, optional
            Filter tags that should exclude the locus from being returned
        filter_requires : set, optional
            Filter tags that should be required for a locus to be returned

        Returns
        -------
        A generator that generates tuples of chrom,pos,{ref,alt1,alt2,...),{field->val}.
        Fields without a corresponding value will have value "None"
        """
        if filter_excludes is None:
            filter_excludes = set()
        if filter_requires is None:
            filter_requires = set()
        # Go to the start of the genotyping part of the file
        file_seek_pos = self.file_handle.seek(self.start_genotype_byte,0)
        # Go through the file
        for line in self.file_handle:
            line = line.lstrip().split('\t',self.header_dict['INFO']+2)
            if len(line) < 2: continue
            # Get chromosome and position
            chrom,pos = map(int,line[:2])
            # Get quality and filter if necessary
            filt = line[self.header_dict['FILTER']].split(',')
            if len(filter_excludes.intersection(filt)) > 0: continue
            elif len(filter_requires.intersection(filt)) != len(filter_requires): continue
            # Get the possible alleles
            alt_alleles = line[self.header_dict['ALT']]
            if alt_alleles != '.':
                alleles = tuple([line[self.header_dict['REF']]]+\
                                line[self.header_dict['ALT']].split(','))
            else:
                alleles = (line[self.header_dict['REF']],)
            ## Create dictionary for info
            info_dict = {}
            for field in line[self.header_dict['INFO']].split(';'):
                if '=' in field:
                    key,val = field.split('=')
                    # Check for regular expressions
                    if int_re.match(val):
                        val = int(val)
                    elif float_re.match(val):
                        val = float(val)
                    elif multi_int_re.match(val):
                        val = tuple(map(int, val.split(',')))
                    elif multi_float_re.match(val):
                        val = tuple((map(float,val.split(','))))
                    info_dict[key] = val
                else:
                    key = field
                    info_dict[key] = None
            yield chrom, pos, alleles, info_dict
            self.current_line += 1
    @staticmethod
    def get_affected_ref_bases(vcf_pos, ref_allele, alt_allele):
        """ Gets the reference positions that are modified through
        the alternative allele, either through substitution or deletion
        
        Parameters
        ----------
        vcf_pos : int
            The position given for this variant in the VCF file
        ref_allele : str
            The reference allele given by the VCF file
        alt_allele : str
            The alternative allele given by the VCF file

        Returns
        -------
        Set of positions that are either substituted or deleted in the
        alternative alleles

        Examples
        --------
        >>> VCF_parser.get_affected_ref_bases(20,'C','T')
        set([20])
        >>> VCF_parser.get_affected_ref_bases(20,'C','CTAG')
        set([])
        >>> VCF_parser.get_affected_ref_bases(20,'TCG','T')
        set([21, 22])
        >>> VCF_parser.get_affected_ref_bases(20,'TCGCG','TCG')
        set([24, 23])
        >>> VCF_parser.get_affected_ref_bases(20,'TCGCG','TCGCGCG')
        set([])
        """
        bases = frozenset(['A','C','G','T'])
        affected_set = set()
        if len(ref_allele) == 1 and len(alt_allele) == 1:
            if alt_allele not in bases:
                raise TypeError("%s is not a valid allele" % alt_allele)
            elif ref_allele != alt_allele:
                return set([vcf_pos])
            else:
                return set()
        elif len(ref_allele) >= len(alt_allele):
            for i in xrange(len(ref_allele)):
                if i < len(alt_allele):
                    # Adding substituted bases
                    if ref_allele[i] != alt_allele[i]:
                        affected_set.add(i+vcf_pos)
                else:
                    # Adding deleted bases
                    affected_set.add(i+vcf_pos)
        elif len(ref_allele) < len(alt_allele):
            for i in xrange(len(ref_allele)):
                # Adding substituted bases
                if ref_allele[i] != alt_allele[i]:
                    affected_set.add(i+vcf_pos)
        return affected_set
            
        
