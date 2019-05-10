
import htsjdk.common.Mutation;
import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author Mario Fruzangohar
 */
public class MUTProcessor {
        String fileName;
        public void buildSampleSeq(int min_flanking, int max_flanking, String ref_file, String refIndex, String genotype_file) throws Exception{
        //written 5/6/2018
        //this method adds a new column, first extract reference and update ref with sample allele that is in genotype_file
        
        Writer writer = null;
        BufferedReader mut = null;
        
        FastaSequenceIndex fsi = null;
        IndexedFastaSequenceFile ifsi = null;
        
            fsi = new FastaSequenceIndex(new File(refIndex));
            ifsi =new IndexedFastaSequenceFile(new File(ref_file) , fsi);
            
            File f = new File(this.fileName);
            String file_name = f.getName().substring(0, f.getName().lastIndexOf('.'));
            String file_ext = f.getName().substring(f.getName().lastIndexOf('.')+1);
            writer = new FileWriter(new File(f.getParent(), file_name+"_added_seq."+file_ext));
            mut = new BufferedReader(new FileReader(this.fileName));
            
            
            String mut_line = mut.readLine();
            while (mut_line != null && !mut_line.isEmpty()){
                String parent_seq = null; 
                Mutation mutation = new Mutation();
                mutation.getFromString8(mut_line);
                int start_idx = 0; 
                int end_idx = 0; 
                int contig_sz = (int)fsi.getContigSize(mutation.contig);
                try{
                    if (mutation.position>=min_flanking && mutation.position<=contig_sz-min_flanking){
                        if (mutation.position>=max_flanking)
                            start_idx = mutation.position - max_flanking;
                        else
                            start_idx = 1;

                        if (mutation.position<=contig_sz-max_flanking)
                            end_idx = mutation.position + max_flanking;
                        else
                            end_idx = contig_sz;
                    }
                    
                    if (start_idx == 0 || end_idx == 0)
                        throw new Exception("Not enough sequences flanking mutation in reference file");
                        ReferenceSequence refseq = ifsi.getSubsequenceAt(mutation.contig, start_idx, end_idx);
                        int mut_idx = mutation.position - start_idx; 
                        parent_seq = buildPolymorphicSeq(genotype_file, mutation.contig, start_idx, refseq.getBaseString(), mutation.position);
                        int asterisk_idx = parent_seq.indexOf('*');
                        if (asterisk_idx<0)
                           parent_seq = parent_seq.substring(0, mut_idx)+mutation.print4website()+parent_seq.substring(mut_idx + 1);
                        else
                            parent_seq = parent_seq.substring(0, asterisk_idx)+mutation.print4website()+parent_seq.substring(asterisk_idx + 2);
                        writer.write(String.format("%s\t%d\t%d\t%s\n", mut_line, start_idx, end_idx,  parent_seq ));
                    }catch(Exception ex){
                        System.out.println(ex.getMessage());
                        writer.write(mut_line + "\n");
                    }
                mut_line = mut.readLine();
            } //while
            
        
        
    }
    private String buildPolymorphicSeq(String genotype_file, String contig, int start_pos, String ref_seq, int mut_pos) throws Exception{
        BufferedReader genotype_reader = null;
        StringBuffer buff = new StringBuffer();
        int end_pos = start_pos + ref_seq.length()-1;
        
            genotype_reader = new BufferedReader(new FileReader(genotype_file));
            String gLine = genotype_reader.readLine();
            String[] arr = gLine.split("\\s");
            String aContig = arr[0].trim();
            int aPos = Integer.parseInt(arr[1]);
            String gStr = arr[3];
            
            while (gLine != null && !aContig.equals(contig)){                
                gLine = genotype_reader.readLine();
                if(gLine == null)  continue;
                arr = gLine.split("\\s");
                aContig = arr[0].trim();
                aPos = Integer.parseInt(arr[1]);
                gStr = arr[3];
            }
            if (gLine == null)
                return ref_seq;
            
            if (aPos > end_pos)
                return ref_seq;
            
            while(gLine != null && aContig.equals(contig) && aPos < start_pos){
                gLine = genotype_reader.readLine();
                if(gLine == null)  continue;
                arr = gLine.split("\\s");
                aContig = arr[0].trim();
                aPos = Integer.parseInt(arr[1]);
                gStr = arr[3];
            }
            
            if(gLine==null || aPos >end_pos || !aContig.equals(contig)) 
                return ref_seq;
            
            int last_pos = start_pos-1;
            while (gLine != null  && aContig.equals(contig) && aPos <= end_pos){
                for (int i = last_pos + 1 ; i < aPos ; i++){
                    if (i==mut_pos)
                        buff.append('*');
                     buff.append(ref_seq.charAt(i-start_pos));
                }
                    
                if (aPos==mut_pos)
                        buff.append('*');

                    int slash_idx = gStr.indexOf('/'); 
                    String allele = null;
                    String insertion = null;
                    if (slash_idx>0){
                       allele = gStr.substring(0, slash_idx);
                       insertion = gStr.substring(slash_idx+1);
                    }else
                        allele = gStr;
                    
                    if(allele.length()==1 && !allele.equals("D")){
                        buff.append(allele);
                    }else if (allele.length()>1){
                        buff.append(String.format("[%s]", allele.replace('D', (char)0)));
                    }
                    if(insertion != null){
                        slash_idx = insertion.indexOf('/');
                        if (slash_idx<0){
                            buff.append(insertion.substring(1)); 
                        }else{
                            String ins1 = insertion.substring(1, slash_idx); 
                            String ins2 = insertion.substring(slash_idx+2); 
                            if (ins1.isEmpty())
                                buff.append(String.format("[%s]", ins2));
                            else if (ins2.isEmpty())
                                buff.append(String.format("[%s]", ins1));
                            else{
                                buff.append(String.format("[%s|%s]", ins1, ins2));
                            }
                        }
                    }
                    
                    
               
                gLine = genotype_reader.readLine();
                if(gLine == null)  continue;
                arr = gLine.split("\\s");
                aContig = arr[0].trim();
                aPos = Integer.parseInt(arr[1]); 
                gStr = arr[3];
            }
            
           
        
        return buff.toString();
    }

}
