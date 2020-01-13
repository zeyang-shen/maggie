import numpy as np
import pandas as pd
import os

import matplotlib.pyplot as plt
plt.switch_backend('agg')
import seaborn as sns


def save_logos(motif_dict, folder='.', input_file='maggie_output_mergedSignificant.tsv'):
    if not os.path.exists(folder+'/logos'):
        os.mkdir(folder+'/logos')
    
    read_df = pd.read_csv(folder+'/'+input_file, sep='\t', index_col=0)
    try:
        import logomaker as lm
    except:
        print('ERROR: missing required package to generate motif logos')
        for k, i in enumerate(read_df.index.values):
            fig = plt.figure(figsize=(6,2))
            fig.savefig(folder+'/logos/'+str(k+1)+'.png', format='png')
            plt.close();
        return
    
    for k, i in enumerate(read_df.index.values):
        mList = i.split('|')
        show_motif = mList[0] # display the motif logo for the top motif in the merged list
        display_df = pd.DataFrame(motif_dict[show_motif].pssm)
        display_df = display_df[['A', 'C', 'G', 'T']]
        display_df.index.name = 'pos'
        display_df = display_df*(display_df>0)

        # plot motif logo
        sns.set(style='white', font_scale=1.5)
        logo = lm.Logo(df=display_df,
                       font_name='DejaVu Sans',
                       fade_below=0.8,
                       shade_below=0.1, figsize=(6,2))
        sns.despine(top=True, bottom=True)
        plt.tight_layout()
        # set axes labels
    #     logo.ax.set_xlabel('Pos')
        logo.ax.set_ylabel('Bits')
        logo.ax.set_title(show_motif)
        logo.ax.set_xticks(np.arange(len(display_df)))
        logo.fig.savefig(folder+'/logos/'+str(k+1)+'.png', format='png')
        plt.close();
    print('Successfully saved motif logos')


def save_distribution(data_df, folder='.', input_file='maggie_output_mergedSignificant.tsv'):
    if not os.path.exists(folder+'/distributions'):
        os.mkdir(folder+'/distributions')
    
    read_df = pd.read_csv(folder+'/'+input_file, sep='\t', index_col=0)
    for k, i in enumerate(read_df.index.values):
        mList = i.split('|')
        show_motif = mList[0].split('$')[1] # plot for the top motif in the merged list
        diffs = np.array(data_df.loc[show_motif, 'score difference'])
        nonzero_diff = diffs[diffs!=0]
        max_val = np.max(np.abs(nonzero_diff))
        # plot distribution
        sns.set(style='whitegrid', font_scale=1.5)
        fig = plt.figure(figsize=(5,2))
        if np.median(nonzero_diff) > 0:
            color = 'salmon'
        else:
            color = 'lightskyblue'
        sns.boxplot(nonzero_diff, color=color)
        plt.xlim(-max_val, max_val)
        plt.axvline(0, c='forestgreen', linestyle='--')
        sns.despine(top=True, bottom=True, left=True)
        plt.tight_layout()
        # set axes labels
        fig.savefig(folder+'/distributions/'+str(k+1)+'.png', format='png')
        plt.close();
    print('Successfully saved distribution plots')


def generate_html(folder='.', input_file='maggie_output_mergedSignificant.tsv'):
    read_df = pd.read_csv(folder+'/'+input_file, sep='\t', index_col=0)
    tot_seq = int(read_df.index.name.split('Total sequences: ')[1][:-1])
    # universal designs
    header = '<head>\n<title>Maggie results</title>\n</head>'
    basic_setup = '<body bgcolor="white" text="black">\n<table border="2" align="center" width="50%" height="20%" bordercolor="grey" cellspacing="5" cellpadding="20">'
    basic_info = 'Total sequences = '+str(tot_seq)
    caption = '<caption><font size="6", color="green"><b>Significant functional motifs</b></font><BR><font size="4", color="black">'+basic_info+'</font></caption>'
    label_color = '#1387FF'
    table_label = '<tr>\n<th width="15%" height="4%"><font color="'+label_color+'">Rank</font></th>\n<th><font color="'+label_color+'">Motif(s)</font></th>\n<th><font color="'+label_color+'">Motif logo</font></th>\n<th>\n<font color="'+label_color+'">Signed -log10 p-value<BR>[90% CI]\n</th>\n<th><font color="'+label_color+'"># mutation (% total seq)</font></th>\n<th><font color="'+label_color+'"># pos mutation (% total mutation)</font></th>\n<th><font color="'+label_color+'"># neg mutation (% total mutation)</font></th>\n<th><font color="'+label_color+'">Median score difference</font></th>\n<th><font color="'+label_color+'">Mean score difference</font></th>\n<th><font color="'+label_color+'">Score difference distribution<BR><i>positive seq. - negative seq.</i></font></th>\n</tr>'
    ending = '</table>\n</body>\n</html>\n'
    with open(folder+'/mergedSignificant.html', 'w') as hf:
        hf.write('<html>\n')
        hf.write('\n'.join([header, basic_setup, caption, table_label]))
        
        # insert rows of significant motifs
        for k, m in enumerate(read_df.index.values):
            # show motif logo
            img_tag = '<img src="logos/'+str(k+1)+'.png" width="350"/>'
            # show score difference distribution
            img_distr = '<img src="distributions/'+str(k+1)+'.png" width="350"/>'
            # show p-value
            pval = str(read_df.loc[m, 'Median p-val'])
            ci = '['+str(read_df.loc[m, '5% p-val'])+', '+str(read_df.loc[m, '95% p-val'])+']'
            # show mutation info
            num_mut = int(read_df.loc[m, 'All mutation'])
            num_mut_block = str(num_mut)+' ('+str(np.around(num_mut/tot_seq*100, decimals=2))+'%)'
            pos_mut = int(read_df.loc[m, 'Pos mutation'])
            neg_mut = int(read_df.loc[m, 'Neg mutation'])
            pos_mut_block = str(pos_mut)+' ('+str(np.around(pos_mut/num_mut*100, decimals=2))+'%)'
            neg_mut_block = str(neg_mut)+' ('+str(np.around(neg_mut/num_mut*100, decimals=2))+'%)'
            # show score difference info
            median_diff = str(np.around(read_df.loc[m, 'Median score difference'], decimals=2))
            mean_diff = str(np.around(read_df.loc[m, 'Mean score difference'], decimals=2))
            # display a row
            one_row = '<tr>\n<th>'+str(k+1)+'</th>\n<th>'+' | '.join([mm.split('$')[0] for mm in m.split('|')])+'</th>\n<th>'+img_tag+'</th>\n<th>\n'+pval+'<BR>\n'+ci+'\n</th>\n<th>'+num_mut_block+'</th>\n<th>'+pos_mut_block+'</th>\n<th>'+neg_mut_block+'</th>\n<th>'+median_diff+'</th>\n<th>'+mean_diff+'</th>\n<th>'+img_distr+'</th>\n</tr>\n'
            hf.write(one_row)
            
        # ending
        hf.write(ending)