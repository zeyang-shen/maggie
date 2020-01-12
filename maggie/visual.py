import numpy as np
import pandas as pd
import os

import matplotlib.pyplot as plt
plt.switch_backend('agg')
import seaborn as sns


def save_logos(motif_dict, folder='.', input_file='maggie_output_mergedSignificant.tsv'):
    try:
        os.mkdir(folder+'/logos')
    except:
        pass
    
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
        logo.ax.set_ylabel('PSSM')
        logo.ax.set_title(show_motif)
        logo.ax.set_xticks(np.arange(len(display_df)))
        logo.fig.savefig(folder+'/logos/'+str(k+1)+'.png', format='png')
        plt.close();
    print('Successfully saved motif logos')


def generate_html(folder='.', input_file='maggie_output_mergedSignificant.tsv'):
    read_df = pd.read_csv(folder+'/'+input_file, sep='\t', index_col=0)
    
    # universal designs
    header = '<head>\n<title>Maggie results</title>\n</head>'
    basic_setup = '<body bgcolor="white" text="black">\n<table border="2" align="center" width="30%" height="40%" bordercolor="grey" cellspacing="5" cellpadding="30">'
    caption = '<caption><font size="6", color="green"><b>Significant functional motifs</b></font></caption>'
    label_color = '#1387FF'
    table_label = '<tr>\n<th width="7%" height="12%"><font color="'+label_color+'">Rank</font></th>\n<th><font color="'+label_color+'">Motif(s)</font></th>\n<th><font color="'+label_color+'">PSSM logo</font></th>\n<th>\n<font color="'+label_color+'">Avg. signed -log10(p)</font><BR>\n<font color="'+label_color+'">[90% CI]</font>\n</th>\n</tr>'
    ending = '</table>\n</body>\n</html>\n'
    with open(folder+'/mergedSignificant.html', 'w') as hf:
        hf.write('<html>\n')
        hf.write('\n'.join([header, basic_setup, caption, table_label]))
        
        # insert rows of significant motifs
        for k, m in enumerate(read_df.index.values):
            img_tag = '<img src="logos/'+str(k+1)+'.png" width="350"/>'
            pval = str(np.round(read_df.loc[m, 'Median p-val'],4))
            ci = '['+str(np.round(read_df.loc[m, '5% p-val'],4))+','+str(np.round(read_df.loc[m, '95% p-val'],4))+']'
            one_row = '<tr>\n<th>'+str(k+1)+'</th>\n<th>'+'|'.join([mm.split('$')[0] for mm in m.split('|')])+'</th>\n<th>'+img_tag+'</th>\n<th>\n'+pval+'<BR>\n'+ci+'\n</th>\n</tr>\n'
            hf.write(one_row)
            
        # ending
        hf.write(ending)