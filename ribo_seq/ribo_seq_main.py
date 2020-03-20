__author__ = 'boris zinshteyn'
"""
Intended for processing of ribosome footprint profiling data from mammalian cells
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42 #leaves most text as actual text in PDFs, not outlines
import os
import sys
import argparse
import subprocess

import ribo_settings
import ribo_utils
import ribo_lib
import ribo_qc
import ribo_plotting
import ribo_tables

class experiment:
    def __init__(self, settings, threads):
        self.threads = threads
        self.settings = settings
        self.all_settings = [lib_setting for lib_setting in self.settings.iter_lib_settings()]
        self.num_datasets = len(self.all_settings)
        self.num_instances = min(self.num_datasets, self.threads)
        self.threads_per_instance = max(self.threads/self.num_instances-1, 1)

        self.genome_num_instances = min(min(self.num_datasets, self.threads), 5)
        self.genome_threads_per_instance = max(self.threads / self.genome_num_instances, 1)

        self.settings.write_to_log('Initializing experiment %s' % self.settings.get_property('experiment_name'))
        self.num_libs = len([x for x in settings.iter_lib_settings()])
        self.make_ncRNA_mapping_index()
        self.settings.write_to_log('loading GTF annotations')
        self.GTF_annotations = ribo_utils.gtf_data(self.settings.get_annotation_GTF_file(), add_3_for_stop=self.settings.get_property('add_3_for_stop'))
        self.settings.write_to_log('loading genome sequence')
        self.genome = ribo_utils.genome_sequence(self.settings.get_genome_sequence_files())
        if self.settings.get_property('transcriptome_mapping_only'):
            #only mapping to the transcriptome, not the genome, so need to generate the genome to map to
            ribo_utils.make_dir(self.settings.get_transcriptome_dir())
            self.GTF_annotations.write_transcript_sequences_to_FASTA(self.settings.get_transcriptome_FASTA(), self.genome)
            self.settings.write_to_log('loading transcriptome sequence')
            self.transcriptome_sequence = ribo_utils.genome_sequence(self.settings.get_transcriptome_FASTA())
        self.make_genome_mapping_index()
        self.deduplicate_reads()
        self.trim_reads()
        self.remove_adaptor()
        self.map_reads_to_ncrna()
        self.map_reads_to_genome()

        self.initialize_libs()
        self.settings.write_to_log('Finished initializing experiment %s\n' % self.settings.get_property('experiment_name'))

    def make_ncRNA_mapping_index(self):
        make_index = False
        if ribo_utils.file_exists(self.settings.get_property('star_ncrna_dir')):
            self.settings.write_to_log('STAR index exists at %s' % self.settings.get_property('star_ncrna_dir'))
            if self.settings.get_property('rebuild_star_index'):
                self.settings.write_to_log('STAR index rebuild forced')
                make_index = True
            else:
                self.settings.write_to_log('using existing STAR index')
        else:
            make_index = True
            ribo_utils.make_dir(self.settings.get_property('star_ncrna_dir'))
        if make_index:
            self.settings.write_to_log('building STAR index')
            command_to_run = 'STAR --runThreadN %d --runMode genomeGenerate --genomeDir %s --genomeFastaFiles %s*.fa --genomeSAsparseD %d  --genomeSAindexNbases %d 1>>%s 2>>%s' % (
                self.threads, self.settings.get_property('star_ncrna_dir'), self.settings.get_ncrna_sequence_dir(), 1, 6,
                self.settings.get_log(), self.settings.get_log())
            self.settings.write_to_log(command_to_run)
            subprocess.Popen(command_to_run, shell=True).wait()
        self.settings.write_to_log('STAR index ready')

    def make_genome_mapping_index(self):
        make_index = False
        if ribo_utils.file_exists(self.settings.get_property('star_genome_dir')):
            self.settings.write_to_log('STAR index exists at %s' % self.settings.get_property('star_genome_dir'))
            if self.settings.get_property('rebuild_star_index'):
                self.settings.write_to_log('STAR index rebuild forced')
                make_index = True
            else:
                self.settings.write_to_log('using existing STAR index')
        else:
            make_index = True
            ribo_utils.make_dir(self.settings.get_property('star_genome_dir'))
        if make_index:
            self.settings.write_to_log('building STAR index')
            if self.settings.get_property('transcriptome_mapping_only'):
                command_to_run = 'STAR --runThreadN %d --runMode genomeGenerate --genomeDir %s --genomeFastaFiles %s --genomeSAsparseD %d 1>>%s 2>>%s' % (
                    self.threads, self.settings.get_property('star_genome_dir'),
                    self.settings.get_transcriptome_FASTA(),
                    self.settings.get_property('star_index_sparsity'),
                    self.settings.get_log(), self.settings.get_log())
            else:
                command_to_run = 'STAR --runThreadN %d --runMode genomeGenerate --genomeDir %s --genomeFastaFiles %s*.fa --genomeSAsparseD %d 1>>%s 2>>%s' % (
                    self.threads, self.settings.get_property('star_genome_dir'), self.settings.get_genome_sequence_dir(),
                    self.settings.get_property('star_index_sparsity'),
                    self.settings.get_log(), self.settings.get_log())
            self.settings.write_to_log(command_to_run)
            subprocess.Popen(command_to_run, shell=True).wait()
        self.settings.write_to_log('STAR index ready')

    def deduplicate_reads(self):
        for lib_settings in self.settings.iter_lib_settings():
            if not lib_settings.deduplicated_reads_exist():
                break
        else: #else clause executes if break did not occur
            self.settings.write_to_log('using existing deduplicated reads')
            return
        self.settings.write_to_log('deduplicating reads with tally')
        ribo_utils.make_dir(self.rdir_path('deduplicated'))
        ribo_utils.parmap(lambda lib_setting: self.deduplicate_reads_one_lib(lib_setting), self.settings.iter_lib_settings(),
                          nprocs=self.num_instances)
        self.settings.write_to_log('done deduplicating reads with tally')

    def deduplicate_reads_one_lib(self, lib_settings):
        lib_settings.write_to_log('deduplicating reads with tally')

        if self.settings.get_property('deduplicate_reads'):
            command = 'tally -i %s -o %s --with-quality 1>>%s 2>>%s' % (lib_settings.get_fastq_gz_file(),
                                                                        lib_settings.get_deduplicated_reads(),
                                                                        lib_settings.get_log(),
                                                                        lib_settings.get_log())
        else:
            #as currently written, this should NEVER trigger, since the deduplicated reads are just redirected to the
            # fastqgz file in ribo_settings.py
            print "something wrong with deduplication settings!!!"
            sys.exit()
            lib_settings.write_to_log('deduplication turned off, just copying input file')
            command = 'cp %s %s 1>>%s 2>>%s' % (lib_settings.get_fastq_gz_file(),
                                                                        lib_settings.get_deduplicated_reads(),
                                                                        lib_settings.get_log(),
                                                                        lib_settings.get_log())
        lib_settings.write_to_log(command)
        subprocess.Popen(command, shell=True).wait()
        lib_settings.write_to_log('deduplicating reads with tally done')

    def trim_reads(self):
        if not self.settings.get_property('force_retrim'):
            for lib_settings in self.settings.iter_lib_settings():
                if not lib_settings.trimmed_reads_exist():
                    break
            else: #else clause executes if break did not occur
                self.settings.write_to_log('using existing trimmed reads')
                return
        else:
            self.settings.write_to_log('read trimming forced')
        self.settings.write_to_log('trimming reads with seqtk')
        ribo_utils.make_dir(self.rdir_path('trimmed'))
        ribo_utils.parmap(lambda lib_setting: self.trim_reads_one_lib(lib_setting, self.threads_per_instance), self.settings.iter_lib_settings(),
                          nprocs=self.num_instances)
        self.settings.write_to_log('done trimming reads with seqtk')

    def trim_reads_one_lib(self, lib_settings, threads_per_instance):
        lib_settings.write_to_log('trimming_reads')
        bases_to_trim = self.settings.get_property('trim_5p')
        command = 'seqtk trimfq -b %d -e 0 %s | pigz -p %d > %s 2>>%s' % (bases_to_trim,
                                                                                 lib_settings.get_deduplicated_reads(),
                                                                                 threads_per_instance,
                                                                                 lib_settings.get_trimmed_reads(),
                                                                                 lib_settings.get_log())
        lib_settings.write_to_log(command)
        subprocess.Popen(command, shell=True).wait()
        lib_settings.write_to_log('trimming_reads done')

    def remove_adaptor(self):
        if not self.settings.get_property('force_retrim'):
            for lib_settings in self.settings.iter_lib_settings():
                if not lib_settings.adaptorless_reads_exist():
                    break
            else:
                self.settings.write_to_log('using existing adaptorless reads')
                return
        else:
            self.settings.write_to_log('adaptor removal forced')
        self.settings.write_to_log('removing adaptors')
        ribo_utils.make_dir(self.rdir_path('adaptor_removed'))
        ribo_utils.parmap(lambda lib_setting: self.remove_adaptor_one_lib(lib_setting, self.threads_per_instance),
                          self.settings.iter_lib_settings(), nprocs=self.num_instances)
        self.settings.write_to_log('done removing adaptors')

    def remove_adaptor_one_lib(self, lib_settings, threads):
        lib_settings.write_to_log('adaptor trimming')
        """
        -x specifies the 3' adaptor to trim from the forward read
        -Q specifies the lowest acceptable mean read quality before trimming
        -l specifies the minimum post-trimming read length
        -L specifies the maximum post-trimming read length
        -o is the output prefix
        --threads specifies number of threads to use
        """
        command_to_run = 'skewer -x %s -Q %d -l %d -L %d -o %s --quiet --mode any --threads %s %s 1>>%s 2>>%s' % (
            self.settings.get_property('adaptor_3p_sequence'),
            self.settings.get_property('sequence_quality_cutoff'),
            self.settings.get_property('min_post_trimming_length'),
            self.settings.get_property('max_post_trimming_length'),
            lib_settings.get_adaptor_trimmed_reads(prefix_only=True),
            threads,
            lib_settings.get_trimmed_reads(),
            lib_settings.get_log(), lib_settings.get_log())
        lib_settings.write_to_log(command_to_run)
        subprocess.Popen(command_to_run, shell=True).wait()
        compression_command = 'pigz -p %d %s-trimmed.fastq' % (threads, lib_settings.get_adaptor_trimmed_reads(prefix_only=True))
        lib_settings.write_to_log(compression_command)
        subprocess.Popen(compression_command, shell=True).wait()
        lib_settings.write_to_log('adaptor trimming done')

    def map_reads_to_ncrna(self):
        """
        map all reads using bowtie
        :return:
        """
        if not self.settings.get_property('force_remapping'):
            for lib_settings in self.settings.iter_lib_settings():
                if not lib_settings.ncrna_mapped_reads_exist():
                    break
            else:
                self.settings.write_to_log('using existing noncoding RNA mapped reads')
                return
        else:
            self.settings.write_to_log('remapping forced')
        self.settings.write_to_log('mapping reads')
        ribo_utils.make_dir(self.rdir_path('ncrna_mapped_reads'))
        ribo_utils.parmap(lambda lib_setting: self.map_one_library_to_ncrna(lib_setting, self.genome_threads_per_instance),
                          [lib_setting for lib_setting in self.settings.iter_lib_settings() if not lib_setting.ncrna_mapping_finished()], nprocs=self.genome_num_instances)
        self.settings.write_to_log( 'finished mapping reads to noncoding RNA')

    def map_one_library_to_ncrna(self, lib_settings, threads):
        lib_settings.write_to_log('mapping_reads')
        command_to_run = 'STAR --runThreadN %d --limitBAMsortRAM 8000000000 --outBAMsortingThreadN %d --genomeDir %s --readFilesIn %s --readFilesCommand gunzip -c ' \
                         '--outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 2 --outSAMmultNmax 1 --outFileNamePrefix %s ' \
                         '--outReadsUnmapped Fastx 1>>%s 2>>%s' %\
                         (threads, threads, self.settings.get_star_ncrna_dir(),
                          lib_settings.get_adaptor_trimmed_reads(),
                          lib_settings.get_ncrna_mapped_reads_prefix(), lib_settings.get_log(), lib_settings.get_log())
        lib_settings.write_to_log(command_to_run)
        subprocess.Popen(command_to_run, shell=True).wait()
        command_to_run = 'samtools index %s' % (lib_settings.get_ncrna_mapped_reads())
        lib_settings.write_to_log(command_to_run)
        subprocess.Popen(command_to_run, shell=True).wait()
        lib_settings.write_to_log('mapping_reads done')

    def map_reads_to_genome(self):
        """
        map all reads using bowtie
        :return:
        """
        if not self.settings.get_property('force_remapping'):
            for lib_settings in self.settings.iter_lib_settings():
                if not lib_settings.genome_mapped_reads_exist():
                    break
            else:
                self.settings.write_to_log('using existing genome-mapped reads')
                return
        else:
            self.settings.write_to_log('remapping forced')
        self.settings.write_to_log('mapping reads to genome')
        ribo_utils.make_dir(self.rdir_path('genome_mapped_reads'))

        ribo_utils.parmap(lambda lib_setting: self.map_one_library_to_genome(lib_setting, self.genome_threads_per_instance),
                          self.settings.iter_lib_settings(), nprocs=self.genome_num_instances)
        self.settings.write_to_log('finished mapping reads to genome')

    def map_one_library_to_genome(self, lib_settings, threads):
        lib_settings.write_to_log('mapping_reads')

        if self.settings.get_property('transcriptome_mapping_only'):
            command_to_run = 'STAR --runThreadN %d --limitBAMsortRAM 80000000000 --outBAMsortingThreadN %d --genomeDir %s --readFilesIn %s ' \
                             '--outSAMtype BAM SortedByCoordinate --alignEndsType %s ' \
                             '--outFilterType BySJout --outFilterMultimapNmax %d, --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0, --outFilterMismatchNmax %d --outWigType wiggle read1_5p --outFileNamePrefix %s' \
                             ' --outReadsUnmapped Fastx 1>>%s 2>>%s' %\
                             (threads, threads, self.settings.get_star_genome_dir(),
                              lib_settings.get_ncrna_unmapped_reads(),
                              self.settings.get_property('alignendstype'),
                              self.settings.get_property('outfiltermultimapnmax'),
                              self.settings.get_property('outfiltermismatchnmax'),
                              lib_settings.get_genome_mapped_reads_prefix(), lib_settings.get_log(), lib_settings.get_log())
        else:
            command_to_run = 'STAR --runThreadN %d --limitBAMsortRAM 80000000000 --outBAMsortingThreadN %d --genomeDir %s --readFilesIn %s ' \
                             '--outSAMtype BAM SortedByCoordinate --quantTranscriptomeBan Singleend --quantMode TranscriptomeSAM --alignEndsType %s ' \
                             '--alignSJDBoverhangMin %d --alignIntronMax %d --sjdbGTFfile %s --alignSJoverhangMin %d ' \
                             '--outFilterType BySJout --outFilterMultimapNmax %d, --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0, --outFilterMismatchNmax %d --outWigType wiggle read1_5p --outFileNamePrefix %s' \
                             ' --outReadsUnmapped Fastx 1>>%s 2>>%s' %\
                             (threads, threads, self.settings.get_star_genome_dir(),
                              lib_settings.get_ncrna_unmapped_reads(),
                              self.settings.get_property('alignendstype'),
                              self.settings.get_property('alignsjdboverhangmin'),
                              self.settings.get_property('alignintronmax'),
                              self.settings.get_property('annotation_gtf_file'),
                              self.settings.get_property('alignsjoverhangmin'),
                              self.settings.get_property('outfiltermultimapnmax'),
                              self.settings.get_property('outfiltermismatchnmax'),
                              lib_settings.get_genome_mapped_reads_prefix(), lib_settings.get_log(), lib_settings.get_log())
        lib_settings.write_to_log(command_to_run)
        subprocess.Popen(command_to_run, shell=True).wait()

        if not self.settings.get_property('transcriptome_mapping_only'):
            # sort transcript-mapped bam file
            command_to_run = 'samtools sort -@ %d --output-fmt BAM -o %s.temp_sorted.bam %s 1>>%s 2>>%s' % (
            threads, lib_settings.get_transcript_mapped_reads(),
            lib_settings.get_transcript_mapped_reads(),
            lib_settings.get_log(), lib_settings.get_log())
            lib_settings.write_to_log(command_to_run)
            subprocess.Popen(command_to_run, shell=True).wait()
            command_to_run = 'mv %s.temp_sorted.bam %s' % (lib_settings.get_transcript_mapped_reads(),
                                                                              lib_settings.get_transcript_mapped_reads())
            lib_settings.write_to_log(command_to_run)
            subprocess.Popen(command_to_run, shell = True).wait()
            command_to_run = 'samtools index %s' % (lib_settings.get_transcript_mapped_reads())
            lib_settings.write_to_log(command_to_run)
            subprocess.Popen(command_to_run, shell = True).wait()
        command_to_run = 'samtools index %s' % (lib_settings.get_genome_mapped_reads())
        lib_settings.write_to_log(command_to_run)
        subprocess.Popen(command_to_run, shell=True).wait()
        lib_settings.write_to_log('mapping_reads done')

    def perform_qc(self):
        self.settings.write_to_log('performing QC')
        qc_engine = ribo_qc.ribo_qc(self, self.settings, self.threads)
        qc_engine.write_trimming_stats_summary()
        qc_engine.plot_rRNA_size_distributions()
        qc_engine.plot_read_annotations_summary()
        qc_engine.plot_read_annotations_summary(mapped_only=True)
        self.settings.write_to_log('done QC')

    def initialize_libs(self):
        self.settings.write_to_log('initializing libraries, counting reads')
        ribo_utils.make_dir(self.rdir_path('transcript_counts'))
        self.libs = []

        ribo_utils.parmap(lambda lib_settings: ribo_lib.assign_tx_reads(self, self.settings, lib_settings), self.settings.iter_lib_settings(), nprocs = self.threads)
        map(lambda lib_settings: self.initialize_lib(lib_settings), self.settings.iter_lib_settings())
        self.settings.write_to_log('initializing libraries, counting reads, done')

    def find_lib_by_sample_name(self, sample_name):
        for lib in self.libs:
            if lib.lib_settings.sample_name == sample_name:
                return lib
        assert False #if this triggers, your settings file is broken.

    def initialize_lib(self, lib_settings):
        lib = ribo_lib.ribo_lib(self.settings, lib_settings)
        self.libs.append(lib)

    def make_tables(self):
        ribo_utils.make_dir(self.rdir_path('tables'))
        ribo_tables.make_readthrough_table(self)
        ribo_tables.make_detailed_readthrough_table(self)
        ribo_tables.transcriptome_features_table(self)
        ribo_tables.make_cds_rpkm_table(self)
        ribo_tables.make_cds_counts_table(self)

    def make_plots(self):

        ribo_utils.make_dir(self.rdir_path('plots'))

        ribo_plotting.plot_start_codon_average(self, up=100, down=100, min_cds_reads=self.settings.get_property('comparison_read_cutoff'),
                                               read_end='3p', read_lengths='all')
        ribo_plotting.plot_stop_codon_average(self, up=100, down=100, min_cds_reads=self.settings.get_property('comparison_read_cutoff'),
                                               read_end='3p', read_lengths='all')

        ribo_plotting.plot_fragment_length_distributions(self)
        ribo_plotting.plot_frame_distributions(self)

        #ribo_plotting.plot_stop_positional_read_lengths(self, up=100, down=100, min_cds_reads=128, read_end='3p')
        #ribo_plotting.plot_start_positional_read_lengths(self, up=100, down=100, min_cds_reads=128, read_end='3p')
        #ribo_plotting.plot_first_exon_positional_read_lengths(self, up=100, down=100, min_cds_reads=128, read_end='3p')

        #ribo_plotting.plot_readthrough_box(self)
        #ribo_plotting.plot_readthrough_box(self, log=True)

    def rdir_path(self, *args):
        return os.path.join(self.settings.get_rdir(), *args)

    def make_counts_table(self, fractional=False):
        """
        write out number of fragments mapping to each TL in each dataset
        :param fractional: if True, replace raw counts with library fraction in reads per million
        :return:
        """
        if fractional:
            summary_file = open(os.path.join(
                self.rdir_path('tables'),
                'rpm.txt'), 'w')
        else:
            summary_file = open(os.path.join(
                self.rdir_path('tables'),
                'raw_counts.txt'), 'w')

        header = 'sequence name\t' + '\t'.join([lib.lib_settings.sample_name for lib in self.libs]) + '\n'
        summary_file.write(header)
        if fractional:
            for sequence_name in self.libs[0].pool_sequence_mappings:
                out_line = '%s\t%s\n' % (sequence_name,
                                         '\t'.join(['%f' % ((10**6)*lib.pool_sequence_mappings[sequence_name].fragment_count/float(lib.total_mapped_fragments)) for lib in self.libs]))
                summary_file.write(out_line)
        else:
            for sequence_name in self.libs[0].pool_sequence_mappings:
                out_line = '%s\t%s\n' % (sequence_name,
                                         '\t'.join(['%f' %
                                                    lib.pool_sequence_mappings[sequence_name].fragment_count
                                                    for lib in self.libs]))
                summary_file.write(out_line)
        summary_file.close()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("settings_file")
    parser.add_argument("--make-tables",
                        help="Makes tables.",
                        action='store_true')
    parser.add_argument("--perform-qc",
                        help="performs quality control analysis.",
                        action='store_true')
    parser.add_argument("--make-plots",
                        help="Makes plots.",
                        action='store_true')
    parser.add_argument("--all-tasks",
                        help="Makes plots, tables",
                        action='store_true')
    parser.add_argument("--threads",
                        help="Max number of processes to use",
                        type = int, default = 8)
    args = parser.parse_args()

    return args

def main():
    """
    """
    args = parse_args()
    settings = ribo_settings.ribo_settings(args.settings_file)
    ribo_experiment = experiment(settings, args.threads)
    print 'experiment ready'

    if args.make_plots or args.all_tasks:
        print 'plots'
        settings.write_to_log('making plots')
        ribo_experiment.make_plots()
        settings.write_to_log('done making plots')

    if args.make_tables or args.all_tasks:
        print 'tables'
        settings.write_to_log('making tables')
        ribo_experiment.make_tables()
        settings.write_to_log('done making tables')

    if args.perform_qc or args.all_tasks:
        print 'performing QC'
        settings.write_to_log('performing QC')
        ribo_experiment.perform_qc()
        settings.write_to_log('done performing QC')

if __name__ == "__main__":
    # execute only if run as a script
    main()