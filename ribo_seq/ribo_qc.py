from collections import defaultdict
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
import os
import ribo_utils
import numpy as np
import pysam
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy

class ribo_qc:
    def __init__(self, experiment, experiment_settings, threads):
        """
        Constructor for Library class
        """
        self.threads = threads
        self.experiment = experiment
        self.experiment_settings = experiment_settings
        self.experiment_settings.write_to_log('initializing QC engine')
        self.get_property = self.experiment_settings.get_property
        self.get_rdir = experiment_settings.get_rdir
        ribo_utils.make_dir(self.experiment.rdir_path('QC'))
        self.genome = self.experiment.genome
        self.full_QC_GTF_annotations = ribo_utils.gtf_data(self.experiment_settings.get_QC_annotation_GTF_file(), add_3_for_stop=self.experiment_settings.get_property('add_3_for_stop'))
        self.lib_QCs = [self.initialize_qc_lib(lib_settings) for lib_settings in self.experiment_settings.iter_lib_settings()]
        #TODO: libraries need not be initiated if the proper intermediate TSV files exist

    def write_trimming_stats_summary(self):
        out_file = open(self.experiment_settings.get_trimming_count_summary(), 'w')
        out_percent_file = open(self.experiment_settings.get_trimming_percent_summary(), 'w')

        header_list = ['trimming stat'] + [single_lib_qc.lib_settings.sample_name for single_lib_qc in self.lib_QCs]
        header_line = '%s\n' % ('\t'.join(header_list))
        out_file.write(header_line)
        out_percent_file.write(header_line)
        for stat_name in ['reads processed','reads filtered out by quality control',
                          'short reads filtered out after trimming by size control',
                          'empty reads filtered out after trimming by size control',
                          'reads available', 'trimmed reads available']:
            if stat_name in single_lib_qc.adaptor_stats.keys():
                out_list = [stat_name] + [str(single_lib_qc.adaptor_stats[stat_name]) for single_lib_qc in self.lib_QCs]
                out_percent_list = [stat_name] +\
                                   [str(100.*single_lib_qc.adaptor_stats[
                                            stat_name]/float(single_lib_qc.adaptor_stats['reads processed']))
                                    for single_lib_qc in self.lib_QCs]
                out_line = '%s\n' % ('\t'.join(out_list))
                out_file.write(out_line)
                out_percent_line = '%s\n' % ('\t'.join(out_percent_list))
                out_percent_file.write(out_percent_line)

        out_file.close()
        out_percent_file.close()

    def plot_rRNA_size_distributions(self):
        dfs = []
        for qc_lib in self.lib_QCs:
            frag_dict = qc_lib.rrna_fragment_lengths
            frag_lengths = sorted(frag_dict.keys())
            frag_length_counts = [frag_dict[length] for length in frag_lengths]
            d = {'fragment length': frag_lengths, '# reads': frag_length_counts,
                 '% reads': 100. * np.array(frag_length_counts) / sum(frag_length_counts),
                 'sample': [qc_lib.lib_settings.sample_name] * len(frag_length_counts)}
            temp_df = pd.DataFrame(data=d)
            dfs.append(temp_df)
        frag_length_df = pd.concat(dfs)
        outname = os.path.join(self.experiment_settings.get_rdir(), 'QC', 'rRNA_fragment_sizes.tsv')
        frag_length_df.to_csv(outname, sep='\t', index_label=False)
        sns.set(style="white", color_codes=True, font_scale=1)
        sns.set_style("ticks", {"xtick.major.size": 3, "xtick.minor.size": 2, "ytick.major.size": 2})
        g = sns.factorplot(x='fragment length', y='% reads', hue="sample", data=frag_length_df, legend_out=True,
                           kind="point", size=8)
        outname = os.path.join(self.experiment_settings.get_rdir(), 'QC', 'rRNA_fragment_sizes_single_plot.pdf')
        g.savefig(outname, transparent=True)
        g = sns.factorplot(x='fragment length', y='% reads', col="sample", col_wrap=3, data=frag_length_df, legend_out=True,
                           kind="point", size=8)
        outname = os.path.join(self.experiment_settings.get_rdir(), 'QC', 'rRNA_fragment_sizes_multi_plot.pdf')
        g.savefig(outname, transparent=True)

    def plot_read_annotations_summary(self, mapped_only=False, representation_cutoff = 1.0):
        dfs = []
        for qc_lib in self.lib_QCs:
            dict_list = [] # a list of tuples that I wil later cast to a dict
            dict_list.append(('sample',qc_lib.lib_settings.sample_name))
            dict_list.append(('input reads', qc_lib.input_reads))
            dict_list.append(('reads processed', qc_lib.adaptor_stats['reads processed']))
            dict_list.append(('<%d nt' % (self.get_property('min_post_trimming_length')), qc_lib.adaptor_stats['reads processed'] - qc_lib.adaptor_stats['trimmed reads available']))
            for ncrna_reference in qc_lib.ncrna_reference_counts:
                dict_list.append((ncrna_reference, qc_lib.ncrna_reference_counts[ncrna_reference]))
            for uniqueness in qc_lib.annotation_mapping_counts:
                for annotation_type in qc_lib.annotation_mapping_counts[uniqueness]:
                    dict_list.append((uniqueness+'_'+annotation_type, qc_lib.annotation_mapping_counts[uniqueness][annotation_type]))
            temp_df = pd.DataFrame(data=dict(dict_list), index=[qc_lib.lib_settings.sample_name])
            dfs.append(temp_df)
        read_summary = pd.concat(dfs)
        outname = os.path.join(self.experiment_settings.get_rdir(), 'QC', 'read_annotation_summary.tsv')
        read_summary.to_csv(outname, sep='\t')

        #now lets assemble the info we'll need for making the stacked bar graph.
        #first consolidate the tRNA data
        consolidated_summary = pd.DataFrame()
        tRNA_headers = []
        rRNA_headers = []
        snRNA_headers = []

        for header in list(read_summary.columns.values):
            if 'tRNA' in header or 'MT' in header or '(' in header:
                tRNA_headers.append(header)
            elif 'rRNA' in header or 'RDN' in header:
                rRNA_headers.append(header)
            elif 'snR' in header:
                snRNA_headers.append(header)
            else:
                consolidated_summary[header] = read_summary[header]

        consolidated_summary['tRNA_mapping'] = read_summary[tRNA_headers].sum(axis=1)
        consolidated_summary['rRNA_mapping'] = read_summary[rRNA_headers].sum(axis=1)
        consolidated_summary['snRNA_mapping'] = read_summary[snRNA_headers].sum(axis=1)

        outname = os.path.join(self.experiment_settings.get_rdir(), 'QC', 'consolidated_read_annotation_summary.tsv')
        consolidated_summary.to_csv(outname, sep='\t')
        consolidated_percent_summary = consolidated_summary[
            [header for header in sorted(list(consolidated_summary.columns.values)) if header != 'sample']].div(
            consolidated_summary['reads processed'], axis='index')
        consolidated_percent_summary = consolidated_percent_summary.apply(lambda x: x * 100.)
        consolidated_percent_summary['sample'] = consolidated_summary['sample']
        outname = os.path.join(self.experiment_settings.get_rdir(), 'QC', 'consolidated_percent_read_annotation_summary.tsv')
        consolidated_percent_summary.to_csv(outname, sep='\t')
        top_summary = pd.DataFrame()
        totals = consolidated_percent_summary.sum()  # sum each column
        multi_totals = {}
        unique_totals = {}
        for header in list(consolidated_percent_summary.columns.values):
            sp = header.split(' mapping')
            if sp[0] == 'multiple':
                multi_totals[header] = totals[header]
            else:
                unique_totals[header] = totals[header]
        #sort into annotations that are represented above a certain percentage
        top_multi = [anno for anno in sorted(multi_totals.keys(), key=lambda x: multi_totals[x], reverse=True) if max(consolidated_percent_summary[anno])>=representation_cutoff]
        other_multi = [anno for anno in sorted(multi_totals.keys(), key=lambda x: multi_totals[x], reverse=True) if max(consolidated_percent_summary[anno])<representation_cutoff]
        top_unique = [anno for anno in sorted(unique_totals.keys(), key=lambda x: unique_totals[x], reverse=True) if max(consolidated_percent_summary[anno])>=representation_cutoff]
        other_unique = [anno for anno in sorted(unique_totals.keys(), key=lambda x: unique_totals[x], reverse=True) if max(consolidated_percent_summary[anno])<representation_cutoff]
        for header in list(consolidated_percent_summary.columns.values):
            if 'mapping' in header or '<' in header:
                if header in top_multi or header in top_unique:
                    top_summary[header] = consolidated_percent_summary[header]
            else:
                top_summary[header] = consolidated_percent_summary[header]

        top_summary['multiple other'] = consolidated_percent_summary[other_multi].sum(axis=1)
        top_summary['unique other'] = consolidated_percent_summary[other_unique].sum(axis=1)

        top_summary['unmapped'] = top_summary['reads processed'] - consolidated_percent_summary[
            [header for header in list(consolidated_percent_summary.columns.values) if header not in ['reads processed', 'input reads']]].sum(
            axis=1)
        del top_summary['reads processed']
        del top_summary['input reads']

        outname = os.path.join(self.experiment_settings.get_rdir(), 'QC', 'top_read_annotation_summary.tsv')
        top_summary.to_csv(outname, sep='\t')

        if mapped_only:
            del top_summary['unmapped']
            if '<%d nt' % (self.get_property('min_post_trimming_length')) in top_summary:
                del top_summary['<%d nt' % (self.get_property('min_post_trimming_length'))]
            top_summary = top_summary[
                [header for header in sorted(list(top_summary.columns.values)) if header != 'sample']].div(
                top_summary.sum(axis=1), axis='index')
            top_summary = top_summary.apply(lambda x: x * 100.)
            top_summary['sample'] = consolidated_percent_summary['sample']

        # now, we will want to sort these by abundance again
        avg = top_summary.mean()
        avg.sort_values(inplace=True, ascending=False)
        ind = numpy.arange(len(top_summary))  # the x locations for the groups
        width = 0.8  # the width of the bars: can also be len(x) sequence
        plotLayers = []
        bottoms = [0] * len(top_summary)
        bottoms = numpy.array(bottoms)

        fig = plt.figure(figsize=(8, 8))
        num_plots_wide = 1
        num_plots_high = 1
        plot = fig.add_subplot(num_plots_high, num_plots_wide, 1)

        hatch_index = 1
        unhatch_index = 1
        for anno_type, val in avg.iteritems():
            if 'multiple' in anno_type:
                plotLayers.append(plot.bar(ind, top_summary[anno_type].values, width, bottom=bottoms,
                                           color=ribo_utils.colors[hatch_index % len(ribo_utils.colors)],
                                           label=anno_type, hatch='////'))
                hatch_index += 1
            else:
                plotLayers.append(plot.bar(ind, top_summary[anno_type].values, width, bottom=bottoms,
                                           color=ribo_utils.colors[unhatch_index % len(ribo_utils.colors)],
                                           label=anno_type, hatch=None))
                unhatch_index += 1
            bottoms = bottoms + numpy.array(top_summary[anno_type].values)

        plt.ylabel('percent of reads')
        plt.title('summary of read loss in pipeline')
        plot.set_xticks(ind)
        plot.set_xticklabels(top_summary['sample'].values, rotation=85)
        plt.tight_layout()
        # Shrink current axis by 40%
        box = plot.get_position()
        plot.set_position([box.x0, box.y0, box.width * 0.6, box.height])
        # Put a legend to the right of the current axis
        plot.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plot.set_ylim(0, 100)
        if mapped_only:
            outName=os.path.join(self.experiment_settings.get_rdir(), 'QC', 'read_loss_summary_mapped.pdf')
        else:
            outName = os.path.join(self.experiment_settings.get_rdir(), 'QC', 'read_loss_summary.pdf')
        # plt.subplots_adjust(bottom=0.38, right=0.8)
        plt.savefig(outName, transparent=True)

    def initialize_qc_lib(self, lib_settings):
        '''
        if lib_settings.qc_pickle_exists():
            lib_settings.write_to_log('existing QC counts found, unpickling %s' % lib_settings.get_qc_pickle())
            return ribo_utils.unPickle(lib_settings.get_qc_pickle())
        else:
            lib_settings.write_to_log('no existing QC counts found, counting')
            lib = single_lib_qc(self, lib_settings)
            lib_settings.write_to_log('writing %s' % (lib_settings.get_qc_pickle()))
            ribo_utils.makePickle(lib, lib_settings.get_qc_pickle())
        '''
        return single_lib_qc(self, lib_settings)

class single_lib_qc():
    def __init__(self, parent_qc, lib_settings):
        """
        Constructor for Library class 
        """
        self.parent_qc = parent_qc
        self.lib_settings = lib_settings
        self.input_reads = ribo_utils.file_len(lib_settings.get_fastq_gz_file())/4
        self.get_adaptor_trimming_stats()
        self.ncrna_samfile = pysam.AlignmentFile(self.lib_settings.get_ncrna_mapped_reads(), "rb")
        self.count_ncrna_mapping_reads()
        self.ncrna_samfile.close()
        self.genome_samfile = pysam.AlignmentFile(self.lib_settings.get_genome_mapped_reads(), "rb")
        self.get_mapping_multiplicity_stats()
        self.get_mapping_annotation_summary(self.genome_samfile)
        self.genome_samfile.close()

    def count_ncrna_mapping_reads(self):
        self.ncrna_reference_counts = defaultdict(int)
        self.ncrna_sequence_multiplicities = {}
        self.rrna_fragment_lengths = defaultdict(int)
        self.total_rRNA_alignments = 0
        self.lib_settings.write_to_log('counting ncrna-mapped reads')
        for reference in self.ncrna_samfile.references:
                for alignment in self.ncrna_samfile.fetch(reference=reference):
                    if not alignment.is_secondary:
                        self.ncrna_reference_counts[reference] += 1
                        if 'rRNA' in reference:
                            self.total_rRNA_alignments += 1
                            self.rrna_fragment_lengths[alignment.query_length] += 1
                        if alignment.query_sequence in self.ncrna_sequence_multiplicities:
                            self.ncrna_sequence_multiplicities[alignment.query_sequence]['count'] += 1
                        else:
                            self.ncrna_sequence_multiplicities[alignment.query_sequence] = {}
                            self.ncrna_sequence_multiplicities[alignment.query_sequence]['reference'] = reference
                            self.ncrna_sequence_multiplicities[alignment.query_sequence]['count'] = 1
        out_file = open(self.lib_settings.get_ncrna_most_common_reads(), 'w')

        for sequence in sorted(self.ncrna_sequence_multiplicities.keys(), key=lambda x:self.ncrna_sequence_multiplicities[x]['count'], reverse=True)[:100]:
            out_file.write('>%s_%d\n%s\n' % (self.ncrna_sequence_multiplicities[sequence]['reference'], self.ncrna_sequence_multiplicities[sequence]['count'], sequence, ))
        out_file.close()
        self.lib_settings.write_to_log('done counting rrna-mapped reads')


    def get_adaptor_trimming_stats(self):
        """
        Parses the outpur from skewer to consolidate stats
        :param lib_settings: 
        :return: 
        """
        self.lib_settings.write_to_log('parsing adaptor trimming stats in %s' % (self.lib_settings.get_adaptor_trimming_log()))
        self.adaptor_stats = {}
        trimming_log = open(self.lib_settings.get_adaptor_trimming_log())
        for line in trimming_log:
            if "reads processed; of these:" in line:
                self.adaptor_stats['reads processed'] = int(line.strip().split()[0])
            elif "reads filtered out by quality control" in line:
                self.adaptor_stats['reads filtered out by quality control'] = int(line.strip().split()[0])
            elif "short reads filtered out after trimming by size control" in line:
                self.adaptor_stats['short reads filtered out after trimming by size control'] = int(line.strip().split()[0])
            elif "empty reads filtered out after trimming by size control" in line:
                self.adaptor_stats['empty reads filtered out after trimming by size control'] = int(line.strip().split()[0])
            elif "reads available; of these:" in line:
                self.adaptor_stats['reads available'] = int(line.strip().split()[0])
            elif "untrimmed reads available after processing" in line:
                self.adaptor_stats['untrimmed reads available'] = int(line.strip().split()[0])
            elif " trimmed reads available after processing" in line:
                self.adaptor_stats['trimmed reads available'] = int(line.strip().split()[0])
            elif "Length distribution of reads after trimming:" in line:
                self.adaptor_stats['post_trimming_lengths'] = {}
                line = trimming_log.next()
                while True:
                    try:
                        line = trimming_log.next()
                        ll = line.strip().split()
                        length = int(ll[0])
                        count = int(ll[1])
                        self.adaptor_stats['post_trimming_lengths'][length] = count
                    except:
                        break
        trimming_log.close()
        self.lib_settings.write_to_log('done parsing adpator trimming stats')

    def get_mapping_multiplicity_stats(self):
        total_alignments = 0
        primary_alignments = 0
        secondary_alignments = 0
        unmapped_alignments = 0
        multiplicity = defaultdict(int)
        for alignment in self.genome_samfile.fetch():
            if alignment.is_secondary:
                secondary_alignments += 1
            else:
                primary_alignments += 1
                multiplicity[int(alignment.get_tag('NH:i'))] += 1
            if alignment.is_unmapped:
                unmapped_alignments += 1
            total_alignments += 1
        print self.lib_settings.sample_name
        print 'total: ', total_alignments
        print 'primary: ', primary_alignments
        print 'unique mapping: ', multiplicity[1]
        multimapping = sum([multiplicity[mult] for mult in multiplicity.keys() if mult !=1])
        print 'multiply mapping: ', multimapping
        print 'secondary: ', secondary_alignments
        print 'unmapped: ', unmapped_alignments


        #for mult in sorted(multiplicity.keys()):
        #    print '%d: %d' % (mult, multiplicity[mult])

    def get_mapping_annotation_summary(self, genome_samfile, reversed_reads=False):
        '''
        produce a summary of where all of the reads in a dataset map.
        :return: 
        '''
        self.lib_settings.write_to_log('making mapping annotation_summary for all mapping reads')
        self.annotation_mapping_counts = defaultdict(lambda : defaultdict(int))
        self.total_primary_alignments = 0
        for chromosome in genome_samfile.references:
            for alignment in genome_samfile.fetch(reference=chromosome):
                if not alignment.is_secondary:
                    self.total_primary_alignments += 1
                    multiplicity = int(alignment.get_tag('NH:i'))
                    assert multiplicity > 0
                    if multiplicity > 1:
                        uniqueness = 'multiple mapping'
                    else:
                        uniqueness = 'unique mapping'
                    if not self.parent_qc.experiment_settings.get_property('reads_reversed'):#for example, Illumina Truseq
                        if not alignment.is_reverse:  # alignment on + strand
                            strand = '+'
                        else:  # alignment on - strand
                            strand = '-'
                    else:
                        if alignment.is_reverse:  # alignment on + strand
                            strand = '+'
                        else:  # alignment on - strand
                            strand = '-'
                    annotation_entry = self.parent_qc.full_QC_GTF_annotations.find_smallest_annotation_at_position(chromosome, strand, alignment.reference_start, alignment.reference_end)
                    if annotation_entry is None:
                        self.annotation_mapping_counts[uniqueness]['not annotated']+=1
                    elif annotation_entry.get_value('type') is None:
                        self.annotation_mapping_counts[uniqueness]['not annotated'] += 1
                    elif annotation_entry.get_value('type') in ['transcript', 'exon']:
                        if annotation_entry.get_value('transcript_type') is None:
                            self.annotation_mapping_counts[uniqueness][annotation_entry.get_value('type')]+=1
                        else:
                            #if annotation_entry.get_value('transcript_type') == 'protein_coding':
                            #    self.lib_settings.write_to_log('protein_coding hit: %s %s %d %d %s %s' % (chromosome, strand, alignment.reference_start, alignment.reference_end, alignment.cigarstring, alignment.query_sequence))
                            #    self.lib_settings.write_to_log(annotation_entry)
                            self.annotation_mapping_counts[uniqueness][annotation_entry.get_value('transcript_type')]+=1
                    elif annotation_entry.get_value('type') == 'UTR':
                        # need to differentiate 5' from 3' UTR
                        type = self.parent_qc.full_QC_GTF_annotations.utr_type(annotation_entry)
                        assert type is not None
                        self.annotation_mapping_counts[uniqueness][type] += 1
                    else:
                        self.annotation_mapping_counts[uniqueness][annotation_entry.get_value('type')]+=1
                        #if annotation_entry.get_value('ype') == 'protein_coding':
                            #self.lib_settings.write_to_log('protein_coding hit: %s %s %d %d %s %s' % (
                            #chromosome, strand, alignment.reference_start, alignment.reference_end, alignment.cigarstring,
                            #alignment.query_sequence))

        self.lib_settings.write_to_log('total_aligned: %d' % self.total_primary_alignments)
        self.lib_settings.write_to_log('finished making mapping annotation_summary for all mapping reads')