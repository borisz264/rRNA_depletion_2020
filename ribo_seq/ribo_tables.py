import ribo_utils
import numpy as np
import os
import scipy.stats as stats
import math

def make_readthrough_table(experiment):
    all_genes = set()
    sample_names = []
    for lib in experiment.libs:
        sample_name = lib.lib_settings.sample_name
        sample_names.append(sample_name)
        tx_w_data = set([tx.tx_id for tx in lib.transcripts.values() if not tx.compute_readthrough_ratio(16, read_end='3p',
                                                                                           read_lengths='all',
                                                                                           cds_cutoff=128) == None ])
        all_genes = all_genes.union(tx_w_data)
    out_name = os.path.join(experiment.settings.get_rdir(), 'tables', 'readthrough_fractions.tsv')
    f = open(out_name, 'w')
    f.write('tx_id\t%s\n' % '\t'.join(sample_names))
    for tx_name in all_genes:
        values = [str(lib.transcripts[tx_name].compute_readthrough_ratio(16, read_end='3p',read_lengths='all', cds_cutoff=128, log = False)) if
                  tx_name in lib.transcripts and lib.transcripts[tx_name].compute_readthrough_ratio(16, read_end='3p',
                                                                                                    read_lengths='all', cds_cutoff=128, log = False) != None
                  else '' for lib in experiment.libs]
        f.write('%s\t%s\n' % (tx_name, '\t'.join(values)))
    f.close()

def make_cds_rpkm_table(experiment):
    all_genes = set()
    sample_names = []
    for lib in experiment.libs:
        sample_name = lib.lib_settings.sample_name
        sample_names.append(sample_name)
        tx_w_data = set([tx.tx_id for tx in lib.transcripts.values() if tx.is_coding])
        all_genes = all_genes.union(tx_w_data)
    out_name = os.path.join(experiment.settings.get_rdir(), 'tables', 'cds_rpkms.tsv')
    f = open(out_name, 'w')
    f.write('tx_id\tgene_name\t%s\n' % '\t'.join(sample_names))
    for tx_name in all_genes:
        gene_name = lib.get_transcript(tx_name).common_name
        try:
            values = [str(lib.get_cds_rpkm(tx_name, 30, 0, read_end='3p', read_lengths='all')) for lib in experiment.libs]
        except:
            print('CDS too short to compute RPKM with given offsets:', tx_name, tx_name, lib.get_transcript(tx_name).cds_length)
        f.write('%s\t%s\t%s\n' % (tx_name, gene_name, '\t'.join(values)))
        #except:
        #    pass
    f.close()

def make_cds_counts_table(experiment):
    all_genes = set()
    sample_names = []
    for lib in experiment.libs:
        sample_name = lib.lib_settings.sample_name
        sample_names.append(sample_name)
        tx_w_data = set([tx.tx_id for tx in lib.transcripts.values() if tx.is_coding])
        all_genes = all_genes.union(tx_w_data)
    out_name = os.path.join(experiment.settings.get_rdir(), 'tables', 'cds_counts.tsv')
    f = open(out_name, 'w')
    f.write('tx_id\tgene_name\t%s\n' % '\t'.join(sample_names))
    for tx_name in all_genes:
        gene_name = lib.get_transcript(tx_name).common_name
        values = [str(lib.get_transcript(tx_name).get_cds_read_count(30, 0, read_end='3p', read_lengths='all')) for lib in experiment.libs]
        f.write('%s\t%s\t%s\n' % (tx_name, gene_name, '\t'.join(values)))
    f.close()

def make_detailed_readthrough_table(experiment, p_site_offset = 16, read_end='5p', read_lengths='all', cds_cutoff=128,
                                    log=False, post_cds_start_buffer=12, pre_cds_stop_buffer=15,
                                    pre_extension_stop_buffer=15, post_cds_stop_buffer=9):
    #a much more detailed table, which has the raw read numbers and region sizes
    all_genes = set()
    sample_names = []
    p_offset = 16
    headers = []
    for lib in experiment.libs:
        sample_name = lib.lib_settings.sample_name
        sample_names.append(sample_name)
        header_items = ['rt_ratio', 'rt_counts', 'cds_counts']
        for item in header_items:
            headers.append('%s_%s' % (sample_name, item))
        tx_w_data = set([tx.tx_id for tx in lib.transcripts.values()
                         if tx.compute_readthrough_ratio(p_site_offset, read_end=read_end, read_lengths=read_lengths,
                                                         cds_cutoff=cds_cutoff, log=log, post_cds_start_buffer=post_cds_start_buffer,
                                                         pre_cds_stop_buffer=pre_cds_stop_buffer, pre_extension_stop_buffer=pre_extension_stop_buffer,
                                                         post_cds_stop_buffer=post_cds_stop_buffer) is not None and tx.is_coding ])
        all_genes = all_genes.union(tx_w_data)
    headers.append('cds_length')
    headers.append('readthrough_length')
    out_name = os.path.join(experiment.settings.get_rdir(), 'tables', 'detailed_readthrough_fractions.tsv')
    f = open(out_name, 'w')
    f.write('tx_id\t%s\n' % '\t'.join(headers))
    for tx_name in all_genes:
        values = []
        for lib in experiment.libs:
            #want the following parameters:
            #1: readthroug ratio
            #2: readthrough counts
            #3: CDS counts
            if tx_name in lib.transcripts:
                tx = lib.get_transcript(tx_name)
                rt_ratio = str(tx.compute_readthrough_ratio(p_site_offset, read_end=read_end, read_lengths=read_lengths,
                                                         cds_cutoff=cds_cutoff, log=log, post_cds_start_buffer=post_cds_start_buffer,
                                                         pre_cds_stop_buffer=pre_cds_stop_buffer, pre_extension_stop_buffer=pre_extension_stop_buffer,
                                                         post_cds_stop_buffer=post_cds_stop_buffer))
                second_stop = tx.second_stop_position()#zero indexed to beggining of stop codon
                rt_counts = str(tx.get_readthrough_counts(p_offset=p_offset, read_end=read_end, read_lengths=read_lengths,
                                                          pre_extension_stop_buffer=pre_extension_stop_buffer,
                                                          post_cds_stop_buffer=post_cds_stop_buffer))
                cds_counts = str(tx.get_cds_read_count(post_cds_start_buffer-p_offset, -1*pre_cds_stop_buffer-p_offset,
                                                       read_end=read_end, read_lengths=read_lengths))
                if rt_ratio == None:
                    rt_ratio = ''
            else:
                rt_ratio = ''
                rt_counts = ''
                cds_countd = ''
            values.append(rt_ratio)
            values.append(rt_counts)
            values.append(cds_counts)
        #common parameters
        #1: CDS length
        #2: readthrough region length
        values.append(str(tx.cds_length))
        values.append(str(second_stop-(tx.cds_end-2)))
        f.write('%s\t%s\n' % (tx_name, '\t'.join(values)))
    f.close()


def transcriptome_features_table(experiment):
    first_lib = experiment.libs[0]
    all_tx = set(first_lib.transcripts.values())
    out_name = os.path.join(experiment.settings.get_rdir(), 'tables', 'transcript_features.tsv')
    f = open(out_name, 'w')
    f.write('tx_id\tstop_codon_context\tsecond_stop_codon\tUTR_length\tTL_length\tCDS_length\t'
            'tx_length\textension_nt_length\tUTR_A_percent\tUTR_T_percent\tUTR_C_percent\tUTR_G_percent\n')
    for tx in all_tx:
        if tx.is_coding:
            values=[tx.stop_codon_context(), tx.second_stop_codon(), str(tx.trailer_length), str(tx.leader_length),
                    str(tx.cds_length), str(tx.tx_length), str(tx.readthrough_extension_length()),
                    str(tx.trailer_monomer_fraction('A')), str(tx.trailer_monomer_fraction('T')),
                    str(tx.trailer_monomer_fraction('C')), str(tx.trailer_monomer_fraction('G'))]
            values = [x if x not in [None, 'None'] else '' for x in values]
            f.write('%s\t%s\n' % (tx.tx_id, '\t'.join(values)))
    f.close()
