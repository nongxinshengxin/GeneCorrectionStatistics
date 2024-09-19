#! /usr/bin/env python
# -*- coding: utf-8 -*-
# author: <GaoCY>
# email: <EMAIL>

from collections import defaultdict
import argparse
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
#import re

class GFFcompare():
    def __init__(self,before_gff,after_gff,mode,mRNA_geneid,CDS_geneid,mRNA_transid,CDS_transid,output_path):
        self.before_gff=before_gff
        self.after_gff=after_gff
        self.mode=mode
        self.mRNA_geneid=mRNA_geneid
        self.CDS_geneid=CDS_geneid
        self.mRNA_transid=mRNA_transid
        self.CDS_transid=CDS_transid
        self.output_path=output_path

    ##get gene information from gff file
    def gene_parse_gff(self,file_path):
        genes = defaultdict(lambda: {'mRNA': set(), 'CDS': set()})
        chromosomes=defaultdict(set)
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                feature_type = parts[2]
                attributes = parts[8]
                if feature_type == 'exon':
                    
                    for attribute in attributes.split(';'):
                        if attribute.startswith(f'{self.mRNA_geneid}='):
                            gene_id = attribute.split('=')[1]
                            break
                    try:
                        gene_id
                        genes[gene_id]['mRNA'].add((parts[3], parts[4]))
                        chromosomes[parts[0]].add(gene_id)
                    except NameError:
                        continue
                elif feature_type == 'CDS':
                    for attribute in attributes.split(';'):
                        if attribute.startswith(f'{self.CDS_geneid}='):
                            gene_id = attribute.split('=')[1]
                            break
                    try:
                        gene_id
                        genes[gene_id]['CDS'].add((parts[3], parts[4]))
                    except NameError:
                        continue
        return chromosomes,genes
    
    ###get transcript information from gff file
    def transcript_parse_gff(self,file_path):
        chromosomes=defaultdict(set)
        transcripts = defaultdict(lambda: defaultdict(lambda: {'mRNA': set(), 'CDS': set()}))
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                feature_type = parts[2]
                attributes = parts[8]
                if feature_type == 'exon':
                    
                    for attribute in attributes.split(';'):
                        if attribute.startswith(f'{self.mRNA_geneid}='):
                            gene_id = attribute.split('=')[1]
                            
                        if attribute.startswith(f'{self.mRNA_transid}='):
                            trans_id = attribute.split('=')[1]
                            #break
                    try:
                        gene_id,trans_id
                        transcripts[gene_id][trans_id]['mRNA'].add((parts[3], parts[4]))
                        chromosomes[parts[0]].add(gene_id)
                    except NameError:
                        continue

                elif feature_type == 'CDS':
                    for attribute in attributes.split(';'):
                        if attribute.startswith(f'{self.CDS_geneid}='):
                            gene_id = attribute.split('=')[1]
                        if attribute.startswith(f'{self.CDS_transid}='):
                            trans_id = attribute.split('=')[1]
                            #break
                    try:
                        gene_id,trans_id
                        transcripts[gene_id][trans_id]['CDS'].add((parts[3], parts[4]))
                        chromosomes[parts[0]].add(gene_id)
                    except NameError:
                        continue

        return chromosomes,transcripts
    
    ###compare gff files
    def compare_gff(self):
        if self.mode=='gene': 
            _,before_genes = self.gene_parse_gff(self.before_gff)
            chromosomes_after,after_genes = self.gene_parse_gff(self.after_gff)     
            added_genes = set(after_genes) - set(before_genes)
            removed_genes = set(before_genes) - set(after_genes)
        
            corrected_mRNA = set()
            corrected_CDS = set()
        
            for gene_id in before_genes:
                if gene_id in after_genes:
                    if before_genes[gene_id]['mRNA'] != after_genes[gene_id]['mRNA']:
                        corrected_mRNA.add(gene_id)
                    if before_genes[gene_id]['CDS'] != after_genes[gene_id]['CDS']:
                        corrected_CDS.add(gene_id)
            all_corrected_genes = added_genes.union(corrected_mRNA).union(corrected_CDS)
            corrected_genes_in_per_chromosome = defaultdict(set)    
            for chrom in chromosomes_after:
                for gene_per in chromosomes_after[chrom]:
                    if gene_per in all_corrected_genes:
                        corrected_genes_in_per_chromosome[chrom].add(gene_per)
            return corrected_genes_in_per_chromosome,chromosomes_after,added_genes,removed_genes,corrected_mRNA,corrected_CDS                                                                                   
        elif self.mode=='transcript':
            _,before_transcripts = self.transcript_parse_gff(self.before_gff)
            chromosomes_after,after_transcripts = self.transcript_parse_gff(self.after_gff)
            added_genes = set(after_transcripts) - set(before_transcripts)
            removed_genes = set(before_transcripts) - set(after_trasnscripts)

            before_trans_ids = set()
            for sub_dict in before_transcripts.values():
                before_trans_ids.update(sub_dict.keys())
            after_trans_ids = set()

            for sub_dict in after_transcripts.values():
                after_trans_ids.update(sub_dict.keys())
            added_transcripts = after_trans_ids - before_trans_ids
            removed_transcripts = before_trans_ids - after_trans_ids

            transcripts_corrected_mRNA = set()
            transcripts_corrected_CDS = set()
            gene_corrected_mRNA = set()
            gene_corrected_CDS = set()


            for gene_id in before_transcripts:
                if gene_id in after_transcripts:
                    for trans_id in before_transcripts[gene_id]:
                        if trans_id in after_transcripts[gene_id]:
                            if before_transcripts[gene_id][trans_id]['mRNA'] != after_transcripts[gene_id][trans_id]['mRNA']:
                                transcripts_corrected_mRNA.add(trans_id)
                                gene_corrected_mRNA.add(gene_id)
                            if before_transcripts[gene_id][trans_id]['CDS'] != after_transcripts[gene_id][trans_id]['CDS']:
                                transcripts_corrected_CDS.add(trans_id)
                                gene_corrected_CDS.add(gene_id)
            

            all_corrected_genes = added_genes.union(gene_corrected_mRNA).union(gene_corrected_CDS)
            corrected_genes_in_per_chromosome = defaultdict(set)
            for chrom in chromosomes_after:
                for gene_per in chromosomes_after[chrom]:
                    if gene_per in all_corrected_genes:
                        corrected_genes_in_per_chromosome[chrom].add(gene_per)
            return corrected_genes_in_per_chromosome,chromosomes_after,added_genes,removed_genes,added_transcripts,removed_transcripts,transcripts_corrected_mRNA,transcripts_corrected_CDS,gene_corrected_mRNA,gene_corrected_CDS


        else:
            print('Error: mode should be gene or transcript')
            return


    ###output file
    def output_file(self):
        if self.mode=='transcript':
            corrected_genes_in_per_chromosome,chromosomes_after,added_genes,removed_genes,added_transcripts,removed_transcripts,transcripts_corrected_mRNA,transcripts_corrected_CDS,gene_corrected_mRNA,gene_corrected_CDS= self.compare_gff()
        elif self.mode=='gene':
            corrected_genes_in_per_chromosome,chromosomes_after,added_genes,removed_genes,corrected_mRNA,corrected_CDS = self.compare_gff()
        else:
            print('Error: mode should be gene or transcript')

        ####plot
        ###chromosomes distribution bar plot
        chromosomes = list(chromosomes_after.keys())
        positive_values = [len(corrected_genes_in_per_chromosome[chrom]) if chrom in corrected_genes_in_per_chromosome else 0 for chrom in chromosomes ]
        fold = max(positive_values) // 100 + 1
        negative_values = [i/y*100*fold for i,y in zip(positive_values,[len(chromosomes_after[chrom]) for chrom in chromosomes])]

        # Plot setup
        fig, ax = plt.subplots(figsize=(10, 6))

        # Upper bar plot
        ax.bar(chromosomes, positive_values, color='orange', label='Corected genes')

        

        # Lower bar plot (mirrored)
        ax.bar(chromosomes, [-v for v in negative_values], color='teal', label='Percentage of gene (%)')

        # Adding labels for positive values
        for i, v in enumerate(positive_values):
            ax.text(i, v + 2, str(v), ha='center', va='bottom')

        # Adding labels for mirrored negative values
        # for i, v in enumerate(negative_values):
        #     ax.text(i, -v - 5, str(v), ha='center', va='top')

        # Customize plot
        ax.set_xlabel('Chromosomes')
        ax.set_ylabel('per. of gene (%)')
        maxylim = max(positive_values)+(100-max(positive_values) % 100)
        ax.set_ylim(-100*fold, max(positive_values)+10)
        ax.set_yticks([-100*fold, -100*fold/2, 0, maxylim/2, maxylim])
        ax.set_yticklabels([100,50, 0, maxylim/2,maxylim])
        ax.legend(loc='upper right')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.axhline(0, color='black', linewidth=0.8)

        # Show plot
        # plt.tight_layout()
        # plt.show()
        plt.savefig(f'{self.output_path}/my_figure.pdf', format='pdf')
        plt.close()


        if self.mode=='transcript':
            ####amended genes and amended transcripts distribution bar plot
            bar_width = 0.35  

            group1=[len(added_genes),len(gene_corrected_mRNA),len(gene_corrected_CDS)]
            group2=[len(added_transcripts),len(transcripts_corrected_mRNA),len(transcripts_corrected_CDS)]
            x = range(len(group1))

            fig, ax = plt.subplots()  
            rects1 = ax.bar([i-bar_width/2 for i in x], group1, bar_width, color='orange',label='amended genes')  
            rects2 = ax.bar([i+bar_width/2 for i in x], group2, bar_width,  color='teal',label='amended transcripts')  
            
 
            #ax.set_xlabel('Categories')  
            ax.set_ylabel('Number of genes/transcripts')  
            ax.set_xticks(x)  
            ax.set_xticklabels(('Added', 'Corrected mRNA', 'Corrected CDS'))  
            ax.legend()  
            
            # per bar labels 
            def autolabel(rects):  
                for rect in rects:  
                    height = rect.get_height()  
                    ax.annotate('{}'.format(height),  
                                xy=(rect.get_x() + rect.get_width() / 2, height),  
                                xytext=(0, 3),  # 3 points vertical offset  
                                textcoords="offset points",  
                                ha='center', va='bottom')  
            
            autolabel(rects1)  
            autolabel(rects2) 

            plt.savefig(f'{self.output_path}/my_figure2.pdf', format='pdf')
            plt.close()
        elif self.mode=='gene':
            bar_width = 0.35 
            group=[len(added_genes),len(corrected_mRNA),len(corrected_CDS)]
            x = range(len(group))

            rect=plt.bar([i for i in x], group, bar_width, color=['orange','green','blue'])

            plt.xticks([i for i in x], ('Added', 'Corrected mRNA', 'Corrected CDS'))
            #plt.xlabel('Categories')
            plt.ylabel('Number of genes')
            #plt.legend()
            # per bar labels
            def autolabel(rects):  
                for rect in rects:  
                    height = rect.get_height()  
                    plt.annotate('{}'.format(height),  
                                xy=(rect.get_x() + rect.get_width() / 2, height),  
                                xytext=(0, 3),  # 3 points vertical offset  
                                textcoords="offset points",  
                                ha='center', va='bottom')  
            
            autolabel(rect)  
            plt.savefig(f'{self.output_path}/my_figure2.pdf', format='pdf')
            plt.close()
        else:
            print('Error: mode should be gene or transcript')


        ##text file
        percent_values = [i/y*100 for i,y in zip(positive_values,[len(chromosomes_after[chrom]) for chrom in chromosomes])]
        with open(f'{self.output_path}/corrected_messages.txt', 'w') as file:
            if self.mode=='transcript':
                file.write(f"新增基因数: {len(added_genes)}\n")
                file.write(f"删除基因数: {len(removed_genes)}\n")
                file.write(f"新增转录本数: {len(added_transcripts)}\n")
                file.write(f"删除转录本数: {len(removed_transcripts)}\n")
                file.write(f"校正了mRNA的转录本数: {len(transcripts_corrected_mRNA)}\n")
                file.write(f"修改了CDS的转录本数: {len(transcripts_corrected_CDS)}\n")
                file.write(f"校正了mRNA的基因数: {len(gene_corrected_mRNA)}\n")
                file.write(f"修改了CDS的基因数: {len(gene_corrected_CDS)}\n")
                all_correced_genes=added_genes.union(gene_corrected_mRNA).union(gene_corrected_CDS)
                all_correced_transcripts=added_transcripts.union(transcripts_corrected_mRNA).union(transcripts_corrected_CDS)
                file.write(f"所有校正的基因:{len(all_correced_genes)}\n")
                file.write(f"所有校正的转录本:{len(all_correced_transcripts)}\n")
                
            elif self.mode=='gene':
                file.write(f"新增基因数: {len(added_genes)}\n")
                file.write(f"删除基因数: {len(removed_genes)}\n")
                file.write(f"校正的CDS数: {len(corrected_CDS)}\n")
                all_correced_genes=added_genes.union(corrected_mRNA).union(corrected_CDS)
                file.write(f"所有矫正的基因:{len(all_correced_genes)}\n")
            file.write(f"校正的基因在每个染色体上的数量和比例:\n")
            for i,value in enumerate(chromosomes):
                file.write(f"{value}: {positive_values[i]}\n")
                file.write(f"{value}:{percent_values[i]:.2f}%\n")


        with open(f'{self.output_path}/all_corrected_genes_list.txt', 'w') as file:
            if self.mode=='transcript':
                all_corrected_genes = added_genes.union(gene_corrected_mRNA).union(gene_corrected_CDS)

                for gene_id in all_corrected_genes:
                    file.write(f"{gene_id}\n")
                
            elif self.mode=='gene':
                all_corrected_genes = added_genes.union(corrected_mRNA).union(corrected_CDS)
                for gene_id in all_corrected_genes:
                    file.write(f"{gene_id}\n")
        
        if self.mode=='transcript':
            with open(f'{self.output_path}/added_genes_list.txt', 'w') as file:
                for gene_id in added_genes:
                    file.write(f"{gene_id}\n")
            with open(f'{self.output_path}/added_transcripts_list.txt', 'w') as file:
                for trans_id in added_transcripts:
                    file.write(f"{trans_id}\n")
            with open(f'{self.output_path}/removed_transcripts_list.txt', 'w') as file:
                for trans_id in removed_transcripts:
                    file.write(f"{trans_id}\n")
            with open(f'{self.output_path}/corrected_mRNA_transcripts_list.txt', 'w') as file:
                for trans_id in transcripts_corrected_mRNA:
                    file.write(f"{trans_id}\n")
            with open(f'{self.output_path}/corrected_CDS_transcripts_list.txt', 'w') as file:
                for trans_id in transcripts_corrected_CDS:
                    file.write(f"{trans_id}\n")
            with open(f'{self.output_path}/corrected_mRNA_genes_list.txt', 'w') as file:
                for gene_id in gene_corrected_mRNA:
                    file.write(f"{gene_id}\n")
            with open(f'{self.output_path}/corrected_CDS_genes_list.txt', 'w') as file:
                for gene_id in gene_corrected_CDS:
                    file.write(f"{gene_id}\n")
        elif self.mode=='gene':
            with open(f'{self.output_path}/added_genes_list.txt', 'w') as file:
                for gene_id in added_genes:
                    file.write(f"{gene_id}\n")
            with open(f'{self.output_path}/corrected_mRNA_genes_list.txt', 'w') as file:
                for gene_id in corrected_mRNA:
                    file.write(f"{gene_id}\n")
            with open(f'{self.output_path}/corrected_CDS_genes_list.txt', 'w') as file:
                for gene_id in corrected_CDS:
                    file.write(f"{gene_id}\n")

        else:

            print('Error: mode should be gene or transcript')
                

        







def main():
    parser = argparse.ArgumentParser(description='compare gff files')
    parser.add_argument('--before_gff', type=str, help='before gff file')
    parser.add_argument('--after_gff', type=str, help='after gff file')
    parser.add_argument('--mode', type=str, default='gene', help='gene or transcript')
    parser.add_argument('--mRNA_geneid', type=str, default='gene_id',help='exon_geneid')
    parser.add_argument('--CDS_geneid', type=str, default='gene_id',help='CDS_geneid')
    parser.add_argument('--mRNA_transid', type=str, default='Parent',help='exon_transid')
    parser.add_argument('--CDS_transid', type=str, default='Parent',help='CDS_transid')
    parser.add_argument('--output_path', type=str, help='output_path')
    args = parser.parse_args()
    gffcompare = GFFcompare(args.before_gff,args.after_gff,args.mode,args.mRNA_geneid,args.CDS_geneid,args.mRNA_transid,args.CDS_transid,args.output_path)
    #gffcompare.compare_gff()
    gffcompare.output_file()

if __name__ == '__main__':
    main()