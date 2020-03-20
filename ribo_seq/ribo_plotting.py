import ribo_utils
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42 #leaves most text as actual text in PDFs, not outlines
import os
import uniform_colormaps
import math
from matplotlib.ticker import AutoMinorLocator

def plot_fragment_length_distributions(experiment):
    dfs = []
    for lib in experiment.libs:
        frag_dict = lib.get_all_CDS_fragment_length_counts()
        frag_lengths = sorted(frag_dict.keys())
        frag_length_counts = [frag_dict[length] for length in frag_lengths]
        d = {'fragment length': frag_lengths, '# reads': frag_length_counts,
             '% reads': 100. * np.array(frag_length_counts) / sum(frag_length_counts),
             'sample': [lib.lib_settings.sample_name] * len(frag_length_counts)}
        temp_df = pd.DataFrame(data=d)
        dfs.append(temp_df)
        print lib.lib_settings.sample_name, 'reads: ', sum(frag_length_counts)
    frag_length_df = pd.concat(dfs)
    out_name = os.path.join(experiment.settings.get_rdir(), 'plots', 'fragment_length_distributions.tsv')
    frag_length_df.to_csv(out_name, sep='\t')
    out_name = os.path.join(experiment.settings.get_rdir(), 'plots', 'fragment_length_percent_pivot.tsv')
    frag_length_df.pivot(columns='sample', values ='% reads', index='fragment length').to_csv(out_name, sep='\t')

    fig = plt.figure(figsize=(8, 5))
    plots = []
    plot = fig.add_subplot(111)
    color_index = 0
    group_df = frag_length_df.groupby(['sample'])
    for lib in experiment.libs:
        sample = lib.lib_settings.sample_name
        df = group_df.get_group(sample)
        df.plot(x='fragment length', y='% reads', ax=plot, color=ribo_utils.rainbow[color_index%len(ribo_utils.rainbow)], linestyle=ribo_utils.line_styles[color_index/len(ribo_utils.rainbow)], lw=2, sharex=True,
                sharey=True, label=sample)
        color_index += 1
    plots.append(plot)

    for plot in plots:
        #major_xticks = range(12, 60, 3)
        #plot.set_xticks(major_xticks, minor=False)
        # plot.set_xticklabels(major_tick_labels)
        plot.set_ylabel('% CDS-mapping reads', fontsize=20)
        plot.set_xlabel('fragment length', fontsize=20)
        # Hide the right and top spines
        plot.spines['right'].set_visible(False)
        plot.spines['top'].set_visible(False)
        #plot.set_xlim(12, 57)
        #plot.set_ylim(0, 65)
        # try:
        #    plot.legend_.remove()
        # except:
        #    pass
    out_name = os.path.join(experiment.settings.get_rdir(), 'plots', 'fragment_length_distributions.pdf')
    plt.savefig(out_name, transparent='True', format='pdf')
    plt.clf()

def plot_readthrough_box(experiment, log = False):
    fig = plt.figure(figsize=(8, 8))
    num_libs = len(experiment.libs)
    num_plots_wide = 1
    num_plots_high = 1
    colormap = uniform_colormaps.viridis
    plot_index = 0
    plot = fig.add_subplot(num_plots_high, num_plots_wide, plot_index + 1)
    data = []
    legends = []
    boxprops = dict(linewidth=2, color=ribo_utils.black)
    for lib in experiment.libs:
        sample_name = lib.lib_settings.sample_name
        readthroughs = [tx.compute_readthrough_ratio(16, read_end='5p', read_lengths='all', cds_cutoff=128,
                                    log=log, post_cds_start_buffer=12, pre_cds_stop_buffer=15,
                                    pre_extension_stop_buffer=15, post_cds_stop_buffer=9) for
                        tx in lib.transcripts.values() if (not tx.compute_readthrough_ratio(16, read_end='5p', read_lengths='all', cds_cutoff=128,
                                    log=log, post_cds_start_buffer=12, pre_cds_stop_buffer=15,
                                    pre_extension_stop_buffer=15, post_cds_stop_buffer=9) == None) and tx.is_coding ]
        data.append(readthroughs)
        legends.append('%s (%d)' % (sample_name, len(readthroughs)))
        # note that all but the last bin exclude the right (larger) edge of the bin. So I add an extra bin.
    plot.boxplot(data, notch=True, boxprops=boxprops, autorange=True)
    plot_index += 1
    #plot.set_xlabel("fragment length", fontsize=8)
    if log:
        plot.set_ylabel("log10 readthrough fraction", fontsize=8)
    else:
        plot.set_ylabel("readthrough fraction", fontsize=8)
    plot.set_xticklabels(legends, rotation=40, ha='right')
    #plot.set_xlim(min_x, max(bins))
    #lg = plt.legend(loc=2, prop={'size': 12}, labelspacing=0.2)
    #lg.draw_frame(False)
    plt.tight_layout()
    if log:
        out_name = os.path.join(experiment.settings.get_rdir(), 'plots', 'log_readthrough_box.pdf')
    else:
        out_name = os.path.join(experiment.settings.get_rdir(), 'plots', 'readthrough_box.pdf')
    plt.savefig(out_name, transparent='True', format='pdf')
    plt.clf()

def plot_frame_distributions(experiment, read_lengths = ['all', [28], [29], [30]], read_ends = ['5p', '3p']):
    num_libs = len(experiment.libs)
    num_plots_wide = len(read_lengths) * len(read_ends)
    num_plots_high = num_libs
    fig = plt.figure(figsize=(num_plots_wide, num_plots_high))
    colormap = uniform_colormaps.viridis
    plot_index = 0
    for lib in experiment.libs:
        for read_end in read_ends:
            for read_length in read_lengths:
                plot = fig.add_subplot(num_plots_high, num_plots_wide, plot_index+1)
                sample_name = lib.lib_settings.sample_name
                frame_counts = np.zeros(3)
                offsets = {'5p':-15, '3p': 15}
                offset = offsets[read_end]
                for transcript in lib.transcripts.values():
                    if transcript.is_coding:
                        frame_counts = frame_counts + transcript.get_read_frame_counts(transcript.cds_start+offset, transcript.cds_end+offset, read_end=read_end, read_lengths=read_length)

                # note that all but the last bin exclude the right (larger) edge of the bin. So I add an extra bin.
                bar_corners = (np.arange(3)+0.25)*.5
                bar_width = 0.5*.5
                bar_centers = bar_corners + bar_width/2.0
                plot.bar(bar_corners, frame_counts/sum(frame_counts), width=bar_width, label=sample_name, lw=0,
                         color=ribo_utils.black)
                plot.set_ylim(0, 1)
                plot.set_xticks(bar_centers)
                plot.set_xticklabels([str(n) for n in range(3)])
                #plot.set_xlabel("fragment length", fontsize=8)
                plot.set_title('%s, %s' % (read_end, str(read_length)), fontsize=8)
                if plot_index % num_plots_wide == 0:
                    plot.set_ylabel(sample_name, fontsize=8)
                plt.setp(plot.get_xticklabels(), fontsize=7)
                plt.setp(plot.get_yticklabels(), fontsize=7)
                #plot.set_xlim(min_x, max(bins))
                plot_index += 1
    plt.subplots_adjust(wspace=0.5, hspace=0.5)
    #lg = plt.legend(loc=2, prop={'size': 12}, labelspacing=0.2)
    #lg.draw_frame(False)
    #plt.tight_layout()
    out_name = os.path.join(experiment.settings.get_rdir(), 'plots', 'frame_ distributions.pdf')
    plt.savefig(out_name, transparent='True', format='pdf')
    plt.clf()

def plot_start_codon_average(experiment, up = 100, down = 500, min_cds_reads = 128, read_end='5p', read_lengths='all'):
    num_libs = len(experiment.libs)
    num_plots_wide = 1
    num_plots_high = num_libs
    fig = plt.figure(figsize=(8, 2*num_libs))
    colormap = uniform_colormaps.viridis
    plot_index = 0
    for lib in experiment.libs:
        plot = fig.add_subplot(num_plots_high, num_plots_wide, plot_index + 1)
        sample_name = lib.lib_settings.sample_name
        normed_count_sum = np.zeros(down+up+1)
        inclusion_sum = np.zeros(down + up + 1)
        for transcript in lib.transcripts.values():
            if transcript.is_coding:
                if read_end == '5p':
                    start_offset = -15
                    stop_offset = -12
                elif read_end == '3p':
                    start_offset = 14
                    stop_offset = 18
                cds_reads = transcript.get_cds_read_count(start_offset, stop_offset, read_end=read_end, read_lengths=read_lengths)
                if cds_reads >= min_cds_reads:
                    tx_count, tx_inclusion = transcript.get_CDS_read_counts_array(transcript.cds_start, -1 * up, down, read_end=read_end,
                                                                                  read_lengths=read_lengths)
                    normed_count_sum += tx_count/(float(cds_reads)/transcript.cds_length)
                    inclusion_sum += tx_inclusion
        nt_positions = np.arange(-1*up, down+1)-0.5
        plot.bar(nt_positions, normed_count_sum/inclusion_sum,
                  color=colormap((plot_index - 1) / float(len(experiment.libs))), lw=0, label=sample_name)
        plot.set_title(sample_name, fontsize=8)
        plot_index += 1
        if plot_index == num_libs:
            plot.set_xlabel("relative to CDS start", fontsize=8)
        plot.set_ylabel("average density\n (read %s end)" % (read_end), fontsize=8)
        plot.set_xlim(-1 * up, down)
        plot.set_ylim(0, 8)
        minorLocator = AutoMinorLocator(10)
        plot.xaxis.set_minor_locator(minorLocator)
        plot.get_xaxis().set_tick_params(which='both', direction='out')
        plot.get_yaxis().set_tick_params(which='both', direction='out')
    #lg = plt.legend(loc=2, prop={'size': 12}, labelspacing=0.2)
    #lg.draw_frame(False)
    plt.tight_layout()
    out_name = os.path.join(experiment.settings.get_rdir(), 'plots', 'start_codon_avg_%s_%s.pdf' %(read_end, str(read_lengths)))
    plt.savefig(out_name, transparent='True', format='pdf')
    plt.clf()

def plot_second_stop_positions(experiment, up = 500, down = 100, min_cds_reads = 128, read_end='5p', read_lengths='all'):
    num_libs = len(experiment.libs)
    num_plots_wide = 1
    num_plots_high = num_libs
    fig = plt.figure(figsize=(8, 2*num_libs))
    colormap = uniform_colormaps.viridis
    plot_index = 0
    for lib in experiment.libs:
        sample_name = lib.lib_settings.sample_name
        plot = fig.add_subplot(num_plots_high, num_plots_wide, plot_index + 1)
        second_stop_positions = []
        for transcript in lib.transcripts.values():
            if transcript.is_coding:
                cds_reads = transcript.get_cds_read_count(-15, 12, read_end=read_end, read_lengths=read_lengths)
                if cds_reads >= min_cds_reads:
                    second_stop = transcript.second_stop_position()
                    if not second_stop == None:
                        second_stop_positions.append(second_stop+3-transcript.cds_end)
        bins = range(-10, 200)
        bins.append(1000000)
        hist, bin_edges = np.histogram(second_stop_positions, bins)
        plot.bar(bin_edges[:-1]-0.5, hist,
                  color=colormap((plot_index - 1) / float(len(experiment.libs))), lw=0, label=sample_name)
        plot.set_title(sample_name, fontsize=8)
        plot_index += 1
        if plot_index == num_libs:
            plot.set_xlabel("second in-frame stop relative to CDS stop", fontsize=8)
        plot.set_ylabel("# genes", fontsize=8)
        plot.set_xlim(-10, 200)
        minorLocator = AutoMinorLocator(10)
        plot.xaxis.set_minor_locator(minorLocator)
        plot.get_xaxis().set_tick_params(which='both', direction='out')
        plot.get_yaxis().set_tick_params(which='both', direction='out')
    #lg = plt.legend(loc=2, prop={'size': 12}, labelspacing=0.2)
    #lg.draw_frame(False)
    plt.tight_layout()
    out_name = os.path.join(experiment.settings.get_rdir(), 'plots', 'second_stop_positions.pdf')
    plt.savefig(out_name, transparent='True', format='pdf')
    plt.clf()

def plot_stop_codon_average(experiment, up = 500, down = 100, min_cds_reads = 128, read_end='5p', read_lengths='all'):
    '''
    :param experiment:
    :param up:
    :param down:
    :param min_cds_reads:
    :param read_end:
    :param read_lengths:
    :return:
    '''
    num_libs = len(experiment.libs)
    num_plots_wide = 1
    num_plots_high = num_libs
    fig = plt.figure(figsize=(8, 2*num_libs))
    colormap = uniform_colormaps.viridis
    plot_index = 0
    for lib in experiment.libs:
        plot = fig.add_subplot(num_plots_high, num_plots_wide, plot_index + 1)
        sample_name = lib.lib_settings.sample_name
        normed_count_sum = np.zeros(down+up+1)
        inclusion_sum = np.zeros(down + up + 1)
        for transcript in lib.transcripts.values():
            if transcript.is_coding:
                if read_end == '5p':
                    start_offset = -15
                    stop_offset = -12
                elif read_end == '3p':
                    start_offset = 14
                    stop_offset = 18
                cds_reads = transcript.get_cds_read_count(start_offset, stop_offset, read_end=read_end, read_lengths=read_lengths)
                if cds_reads >= min_cds_reads:
                    tx_count, tx_inclusion = transcript.get_CDS_read_counts_array(transcript.cds_end, -1 * up, down, read_end=read_end,
                                                                                  read_lengths=read_lengths)
                    normed_count_sum += tx_count/(float(cds_reads)/transcript.cds_length)
                    inclusion_sum += tx_inclusion
        nt_positions = np.arange(-1*up, down+1)-0.5
        plot.bar(nt_positions, normed_count_sum/inclusion_sum,
                  color=colormap((plot_index - 1) / float(len(experiment.libs))), lw=0, label=sample_name)
        plot.set_title(sample_name, fontsize=8)
        plot_index += 1
        if plot_index == num_libs:
            plot.set_xlabel("nt relative to stop codon", fontsize=8)
        plot.set_ylabel("average density\n (read %s end)" % (read_end), fontsize=8)
        plot.set_xlim(-1*up, down)
        minorLocator = AutoMinorLocator(10)
        plot.xaxis.set_minor_locator(minorLocator)
        plot.get_xaxis().set_tick_params(which='both', direction='out')
        plot.get_yaxis().set_tick_params(which='both', direction='out')
    #lg = plt.legend(loc=2, prop={'size': 12}, labelspacing=0.2)
    #lg.draw_frame(False)
    plt.tight_layout()
    out_name = os.path.join(experiment.settings.get_rdir(), 'plots', 'stop_codon_avg_%s_%s.pdf' %(read_end, str(read_lengths)))
    plt.savefig(out_name, transparent='True', format='pdf')
    plt.clf()

def plot_first_exon_average(experiment, up = 500, down = 100, min_cds_reads = 128, read_end='5p', read_lengths='all'):
    num_libs = len(experiment.libs)
    num_plots_wide = 1
    num_plots_high = num_libs
    fig = plt.figure(figsize=(8, 2*num_libs))
    colormap = uniform_colormaps.viridis
    plot_index = 0
    for lib in experiment.libs:
        plot = fig.add_subplot(num_plots_high, num_plots_wide, plot_index + 1)
        sample_name = lib.lib_settings.sample_name
        normed_count_sum = np.zeros(down+up+1)
        inclusion_sum = np.zeros(down + up + 1)
        for transcript in lib.transcripts.values():
            if transcript.is_coding:
                if not transcript.get_first_jxn_in_CDS() == None:
                    if read_end == '5p':
                        start_offset = -15
                        stop_offset = -12
                    elif read_end == '3p':
                        start_offset = 14
                        stop_offset = 18
                    cds_reads = transcript.get_cds_read_count(start_offset, stop_offset, read_end=read_end, read_lengths=read_lengths)
                    if cds_reads >= min_cds_reads:
                        tx_count, tx_inclusion = transcript.get_CDS_read_counts_array(transcript.get_first_jxn_in_CDS(), -1 * up, down, read_end=read_end,
                                                                                      read_lengths=read_lengths)
                        normed_count_sum += tx_count/(float(cds_reads)/transcript.cds_length)
                        inclusion_sum += tx_inclusion
        nt_positions = np.arange(-1*up, down+1)-0.5
        plot.bar(nt_positions, normed_count_sum/inclusion_sum,
                  color=colormap((plot_index - 1) / float(len(experiment.libs))), lw=0, label=sample_name)
        plot.set_title(sample_name, fontsize=8)
        plot_index += 1
        if plot_index == num_libs:
            plot.set_xlabel("relative to second exon start", fontsize=8)
        plot.set_ylabel("average density\n (read %s end)" % (read_end), fontsize=8)
        plot.set_xlim(-1*up, down)
        minorLocator = AutoMinorLocator(10)
        plot.xaxis.set_minor_locator(minorLocator)
        plot.get_xaxis().set_tick_params(which='both', direction='out')
        plot.get_yaxis().set_tick_params(which='both', direction='out')
    #lg = plt.legend(loc=2, prop={'size': 12}, labelspacing=0.2)
    #lg.draw_frame(False)
    plt.tight_layout()
    out_name = os.path.join(experiment.settings.get_rdir(), 'plots', 'first_ej_avg_%s_%s.pdf' %(read_end, str(read_lengths)))
    plt.savefig(out_name, transparent='True', format='pdf')
    plt.clf()

def plot_stop_positional_read_lengths(experiment, up = 500, down = 100, min_cds_reads = 128, read_end='5p', read_lengths='all'):
    num_libs = len(experiment.libs)
    num_plots_wide = 1
    num_plots_high = num_libs
    fig = plt.figure(figsize=(8, 2*num_libs))
    colormap = uniform_colormaps.viridis
    plot_index = 0
    for lib in experiment.libs:
        plot = fig.add_subplot(num_plots_high, num_plots_wide, plot_index + 1)
        sample_name = lib.lib_settings.sample_name
        length_sum = np.zeros(down+up+1)
        count_sum = np.zeros(down + up + 1)
        for transcript in lib.transcripts.values():
            if transcript.is_coding:
                if read_end == '5p':
                    start_offset = -15
                    stop_offset = -12
                elif read_end == '3p':
                    start_offset = 14
                    stop_offset = 18
                cds_reads = transcript.get_cds_read_count(start_offset, stop_offset, read_end=read_end, read_lengths=read_lengths)
                if cds_reads >= min_cds_reads:
                    length_sum_array, counts_array = transcript.get_avg_read_lengths_array(transcript.cds_end, -1*up, down,
                                                                                           read_end=read_end)
                    length_sum += length_sum_array
                    count_sum += counts_array
        nt_positions = np.arange(-1*up, down+1)
        plot.plot(nt_positions, length_sum/count_sum,
                  color=colormap((plot_index - 1) / float(len(experiment.libs))), lw=1, label=sample_name)
        plot.set_title(sample_name, fontsize=8)
        plot_index += 1
        if plot_index == num_libs:
            plot.set_xlabel("relative to CDS stop", fontsize=8)
        plot.set_ylabel("avg read length\n (read %s end)" % (read_end), fontsize=8)
        plot.set_xlim(-1*up, down)
        #plot.set_ylim(25, 35)
        minorLocator = AutoMinorLocator(10)
        plot.xaxis.set_minor_locator(minorLocator)
        plot.get_xaxis().set_tick_params(which='both', direction='out')
        plot.get_yaxis().set_tick_params(which='both', direction='out')
    #lg = plt.legend(loc=2, prop={'size': 12}, labelspacing=0.2)
    #lg.draw_frame(False)
    plt.tight_layout()
    out_name = os.path.join(experiment.settings.get_rdir(), 'plots', 'stop_lengths_%s_%s.pdf' %(read_end, str(read_lengths)))
    plt.savefig(out_name, transparent='True', format='pdf')
    plt.clf()

def plot_start_positional_read_lengths(experiment, up = 500, down = 100, min_cds_reads = 128, read_end='5p', read_lengths='all'):
    num_libs = len(experiment.libs)
    num_plots_wide = 1
    num_plots_high = num_libs
    fig = plt.figure(figsize=(8, 2*num_libs))
    colormap = uniform_colormaps.viridis
    plot_index = 0
    for lib in experiment.libs:
        plot = fig.add_subplot(num_plots_high, num_plots_wide, plot_index + 1)
        sample_name = lib.lib_settings.sample_name
        length_sum = np.zeros(down+up+1)
        count_sum = np.zeros(down + up + 1)
        for transcript in lib.transcripts.values():
            if transcript.is_coding:
                if read_end == '5p':
                    start_offset = -15
                    stop_offset = -12
                elif read_end == '3p':
                    start_offset = 14
                    stop_offset = 18
                cds_reads = transcript.get_cds_read_count(start_offset, stop_offset, read_end=read_end, read_lengths=read_lengths)
                if cds_reads >= min_cds_reads:
                    length_sum_array, counts_array = transcript.get_avg_read_lengths_array(transcript.cds_start, -1*up, down,
                                                                                           read_end=read_end)
                    length_sum += length_sum_array
                    count_sum += counts_array
        nt_positions = np.arange(-1*up, down+1)
        plot.plot(nt_positions, length_sum/count_sum,
                  color=colormap((plot_index - 1) / float(len(experiment.libs))), lw=1, label=sample_name)
        plot.set_title(sample_name, fontsize=8)
        plot_index += 1
        if plot_index == num_libs:
            plot.set_xlabel("relative to CDS stop", fontsize=8)
        plot.set_ylabel("avg read length\n (read %s end)" % (read_end), fontsize=8)
        plot.set_xlim(-1*up, down)
        #plot.set_ylim(25, 35)
        minorLocator = AutoMinorLocator(10)
        plot.xaxis.set_minor_locator(minorLocator)
        plot.get_xaxis().set_tick_params(which='both', direction='out')
        plot.get_yaxis().set_tick_params(which='both', direction='out')
    #lg = plt.legend(loc=2, prop={'size': 12}, labelspacing=0.2)
    #lg.draw_frame(False)
    plt.tight_layout()
    out_name = os.path.join(experiment.settings.get_rdir(), 'plots', 'start_lengths_%s_%s.pdf' %(read_end, str(read_lengths)))
    plt.savefig(out_name, transparent='True', format='pdf')
    plt.clf()

def plot_first_exon_positional_read_lengths(experiment, up = 500, down = 100, min_cds_reads = 128, read_end='5p', read_lengths='all'):
    num_libs = len(experiment.libs)
    num_plots_wide = 1
    num_plots_high = num_libs
    fig = plt.figure(figsize=(8, 2*num_libs))
    colormap = uniform_colormaps.viridis
    plot_index = 0
    for lib in experiment.libs:
        plot = fig.add_subplot(num_plots_high, num_plots_wide, plot_index + 1)
        sample_name = lib.lib_settings.sample_name
        length_sum = np.zeros(down+up+1)
        count_sum = np.zeros(down + up + 1)
        for transcript in lib.transcripts.values():
            if transcript.is_coding and not transcript.get_first_jxn_in_CDS() == None:
                if read_end == '5p':
                    start_offset = -15
                    stop_offset = -12
                elif read_end == '3p':
                    start_offset = 14
                    stop_offset = 18
                cds_reads = transcript.get_cds_read_count(start_offset, stop_offset, read_end=read_end, read_lengths=read_lengths)
                if cds_reads >= min_cds_reads:
                    length_sum_array, counts_array = transcript.get_avg_read_lengths_array(transcript.get_first_jxn_in_CDS(), -1*up, down,
                                                                                           read_end=read_end)
                    length_sum += length_sum_array
                    count_sum += counts_array
            nt_positions = np.arange(-1*up, down+1)
            plot.plot(nt_positions, length_sum/count_sum, color=colormap((plot_index - 1) / float(len(experiment.libs))), lw=1, label=sample_name)
            plot.set_title(sample_name, fontsize=8)
        plot_index += 1
        if plot_index == num_libs:
            plot.set_xlabel("relative to CDS stop", fontsize=8)
        plot.set_ylabel("avg read length\n (read %s end)" % (read_end), fontsize=8)
        plot.set_xlim(-1*up, down)
        #plot.set_ylim(25, 35)
        minorLocator = AutoMinorLocator(10)
        plot.xaxis.set_minor_locator(minorLocator)
        plot.get_xaxis().set_tick_params(which='both', direction='out')
        plot.get_yaxis().set_tick_params(which='both', direction='out')
    #lg = plt.legend(loc=2, prop={'size': 12}, labelspacing=0.2)
    #lg.draw_frame(False)
    plt.tight_layout()
    out_name = os.path.join(experiment.settings.get_rdir(), 'plots', 'first_ej_lengths_%s_%s.pdf' %(read_end, str(read_lengths)))
    plt.savefig(out_name, transparent='True', format='pdf')
    plt.clf()

def codon_metaplots(library):
    pass
