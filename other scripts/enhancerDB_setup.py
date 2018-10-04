
# coding: utf-8

import pandas as pd

df = pd.read_csv("/home/m.truglio/PycharmProjects/canvas/resources/enhancers_pre_files/Enhancers_celltypes.txt", sep='\t', header=None, low_memory=False) #The celltypes file sent by Zhen
df.columns = ['EnhID', 'chr', 'start', 'end', 'CellType', 'Source']


df_genes = pd.read_csv("/home/m.truglio/PycharmProjects/canvas/resources/enhancers_pre_files/DataDownload-EnhancerTagretGene.txt", sep='\t', low_memory=False) #The target genes file

df['TargetGene'] = df_genes['TargetGene']

df_conserv = pd.read_csv("/home/m.truglio/PycharmProjects/canvas/resources/enhancers_pre_files/DataDownload_ConversationScore.txt", sep='\t', low_memory=False) #The conservation file

df['ConservationMean'] = df_conserv['Mean']
df['ConservationMedian'] = df_conserv['Median']

out = open('/home/m.truglio/PycharmProjects/canvas/resources/enhancers_pre_files/Fantom5_celltypes_reformat.bed','w')
out.write('chr\tstart\tend\tCellType\n')
with open('/home/m.truglio/PycharmProjects/canvas/resources/enhancers_pre_files/Fantom5_celltypes.bed') as f:
    for line in f:
        if line.startswith("track name"):
            track=line.split('="')[1].split('_differentially')[0]
        else:
            out.write(line.split()[0].split('chr')[1]+'\t'+line.split()[1]+'\t'+line.split()[2]+'\t'+track+'\n')

out.close()
df_fantom = pd.read_csv("/home/m.truglio/PycharmProjects/canvas/resources/enhancers_pre_files/Fantom5_celltypes_reformat.bed", sep='\t', low_memory=False) #The conservation file

print("working")
for index, row in df.iterrows():
    print(index, "out of 2793317")
    if row['Source'] == "FANTOM5":
        # print("This is the row")
        # print(row)
        # print(index)
        # input()
        match = (df_fantom[(df_fantom.start==row['start']) & (df_fantom.end==row['end']) & (df_fantom.chr==row['chr'])])
        # print(match)
        if match.empty:
            # print("EMPTY!")
            row['CellType']="NA"
        else:
            
            # print(match['CellType'])
            # print(match['CellType'].tolist())
            row['CellType'] = ','.join(match['CellType'].tolist())
        df.iloc[index] = row
            
df['chr'] = 'chr'+df['chr']
df.to_csv('/home/m.truglio/PycharmProjects/canvas/resources/BED/enhancers_full.txt', sep='\t', header=False, index=False)


