
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.List;
import htsjdk.common.*;
import htsjdk.samtools.reference.FastaSequenceIndex;
import java.util.ArrayList;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author Mario Fruzangohar
 */
public class PileupProcessor {
    public void callMutation2(String pileupFile, String destDir, int sample_num, int min_number, double minor_ratio, int mutation_threshold) throws Exception{
        
        Writer mutationFileUnique = null;
        Writer mutationFileNonUnique = null;
        PileupReader pr = null;
        htsjdk.common.Pileup pu = null;
        int nonUniqueCount = 0;
       
       
            mutationFileUnique = new FileWriter(new File(destDir, "unique_mutation.mut")  , false);  
            mutationFileNonUnique = new FileWriter(new File(destDir, "non_unique_mutation.mut")  , false);  
            long pbIndex = 0;
            double lastpercent = 0.0;
            long fileSz = new File(pileupFile).length();
            if (pileupFile.endsWith(".gz"))
               pr = new PileupReader(pileupFile, sample_num, (byte)1, false);
            else
                pr = new PileupReader(pileupFile, sample_num, false);
            
            
            pu = pr.getNext();
            while (pu != null){
                pbIndex += pu.srcLine.getBytes().length; // we add number of bytes read
                double newpercent = ((double)pbIndex/fileSz)*100;
                if (newpercent-lastpercent>=1){
                    lastpercent = newpercent;
                    System.out.println("progress %"+(int)lastpercent);
                }

                String[] genotypes = pu.callGenotype2(sample_num, min_number, minor_ratio, false, false); // not using mapq, not add coverage
                
                //now see how many samples have genotype at this reference position
                int nonEmptySamples = 0;
                for(String gn: genotypes){
                    if (gn!=null && !gn.isEmpty())
                        nonEmptySamples++;
                }
              
                List<Mutation> mutations = null;
               
                mutations = pu.callMutationsPerSample2(sample_num-1, min_number, minor_ratio); 
                if(mutations.size() == 1){
                    mutationFileUnique.write(String.format("%s\t%d\t%s\t%d\t%s\t%s\t%d\t%d\n",pu.chr_id, pu.pos, Helper.Base_String(pu.ref_nt), mutations.get(0).sample_id , mutations.get(0).print2(), Helper.ZG_Short_String(mutations.get(0).zygosity) ,  mutations.get(0).read_depth, nonEmptySamples));
                    
                }else if (mutations.size()>1){
                     for (Mutation mu : mutations){
                         nonUniqueCount++;
                         mutationFileNonUnique.write(String.format("%s\t%d\t%s\t%d\t%s\t%s\t%d\t%d\n",pu.chr_id, pu.pos, Helper.Base_String(pu.ref_nt), mu.sample_id , mu.print2(), Helper.ZG_Short_String(mu.zygosity) , mu.read_depth, nonEmptySamples));
                     }
                         
                }
                pu = pr.getNext();               
            }//while
            System.out.println("Number of non-unique Mutations:" + nonUniqueCount);
        
        
    }
    
        public void callGenotype2(String pileupFile, String destDir,int sample_num, int min_number, double minor_ratio, boolean useMAPQ, boolean allPositions) throws Exception{
        //output is enotype file and each line looks like: contig   pos    ref_base   genotype_of_sample
        Writer genotypeFile = null;
        PileupReader pr = null;  
        int total_rows = 0; // total rows in pileup;
        int[] polymorphic_counts= new int[sample_num];
        int[] polymorphic_contain_ref_counts= new int[sample_num];
        int[] ins_counts= new int[sample_num];
        htsjdk.common.Pileup pu = null;
       
            genotypeFile = new FileWriter(new File(destDir)  , false);  
            long pbIndex = 0;
            double lastpercent = 0.0;
            long fileSz = new File(pileupFile).length();

            
            
            if (pileupFile.endsWith(".gz"))
              pr = new PileupReader(pileupFile, sample_num, (byte)1, useMAPQ); 
            else
                pr = new PileupReader(pileupFile, sample_num, useMAPQ);  

            pu = pr.getNext();
           
            
            while (pu != null){
                total_rows++;
                pbIndex += pu.srcLine.getBytes().length;
                double newpercent = ((double)pbIndex/fileSz)*100;
                if (newpercent-lastpercent>=1){
                    lastpercent = newpercent;
                    System.out.println("progress %"+(int)lastpercent);
                }
                
                String[] genotypes = pu.callGenotype2(sample_num, min_number, minor_ratio, useMAPQ, false);
                boolean isPolymorphic = false;
                if (allPositions){
                    StringBuffer buff = new StringBuffer();
                    for (int i=0; i<sample_num; i++){
                        buff.append((genotypes[i]==null?"":genotypes[i]) + "\t"); 
                        if (genotypes[i]!=null){
                            int plus_index = genotypes[i].indexOf('+');
                            if (plus_index>0){
                                ins_counts[i]++;
                                if( genotypes[i].indexOf('+', plus_index+1)>0)
                                    isPolymorphic = true;
                            }
                            String part1 = genotypes[i].split("/")[0];
                            if(part1.length()>1){
                                polymorphic_counts[i]++;
                                isPolymorphic = true;
                                if (part1.indexOf(Helper.Base_String( pu.ref_nt))>=0)
                                    polymorphic_contain_ref_counts[i]++;
                            }
                        }
                    }
                    // format : chr_id, start_pos, reference genotype, variation string, genotype of all samples
                    genotypeFile.write(String.format("%s\t%d\t%s\t%s\n",pu.chr_id, pu.pos, Helper.Base_String(pu.ref_nt), buff.toString())  );                        
                }else if(Helper.isVariationFromReference(genotypes, Helper.Base_String( pu.ref_nt))){                    
                    StringBuffer buff = new StringBuffer();
                    for (int i=0; i<sample_num; i++){
                        buff.append((genotypes[i]==null?"":genotypes[i]) + "\t");
                        if (genotypes[i]!=null){
                            int plus_index = genotypes[i].indexOf('+');
                            if (plus_index>0){
                                ins_counts[i]++;
                                if( genotypes[i].indexOf('+', plus_index+1)>0) 
                                    isPolymorphic = true;
                            }
                            String part1 = genotypes[i].split("/")[0];
                            if(part1.length()>1){
                                polymorphic_counts[i]++;
                                isPolymorphic = true;
                                if (part1.indexOf(htsjdk.common.Helper.Base_String( pu.ref_nt))>=0)
                                    polymorphic_contain_ref_counts[i]++;
                            }
                        }
                    }
                    // format : chr_id, start_pos, reference Allele, variation string, genotype of all samples
                    genotypeFile.write(String.format("%s\t%d\t%s\t%s\n",pu.chr_id, pu.pos, Helper.Base_String( pu.ref_nt), buff.toString())  );
                }
                
                
                pu = pr.getNext();
            }//while
            //now print statistics for each sample
            System.out.println("Total Rows:"+total_rows);
            System.out.println("Insertions:");
            for(int i=0; i<sample_num; i++)
               System.out.print(ins_counts[i]+"\t") ;
            System.out.println();
            System.out.println("Polymorphic Positions:");
            for(int i=0; i<sample_num; i++)
               System.out.print(polymorphic_counts[i]+"\t") ;
            System.out.println();
            System.out.println("Polymorphic Positions that contain reference allele:");
            for(int i=0; i<sample_num; i++)
               System.out.print(polymorphic_contain_ref_counts[i]+"\t") ;
            System.out.println();
            
       
        
    }

     public void disocoverHighCoverage(String pileupFile, String refIndexFile, int min_depth, int tail_length, int merge_gap, String destDir) throws Exception{
        
        Writer writer = null;       
        PileupReader pr = null;
        FastaSequenceIndex fsi = null;
        List<SequenceRegion>  regions = new ArrayList();  
        List<SequenceRegion>  merged_regions = new ArrayList();
        
            fsi = new FastaSequenceIndex(new File(refIndexFile));
            writer = new FileWriter(new File(destDir)  , false);  // bed file
             
            long pbIndex = 0;
            double lastpercent = 0.0;
            long fileSz = new File(pileupFile).length();
            
            if (pileupFile.endsWith(".gz"))
               pr = new PileupReader(pileupFile, 1, (byte)1, false);
            else
                pr = new PileupReader(pileupFile, 1, false);

           while (true){
               long start_pos = 0;
               String start_contig = null;
               long end_pos = 0; 

                htsjdk.common.Pileup pu = pr.getNext();
                while (pu != null ){
                    pbIndex += pu.srcLine.getBytes().length; // we add number of bytes read
                    double newpercent = ((double)pbIndex/fileSz)*100;
                    if (newpercent-lastpercent>=1){
                        lastpercent = newpercent;
                        System.out.println("progress %"+(int)lastpercent);
                    }

                    int depth = pu.triplets.get(0).depth;
                    if (depth>=min_depth)
                        break;

                    pu = pr.getNext();

                }//end while
               
                if (pu!=null){
                   
                   start_pos = pu.pos;
                  
                   start_contig = pu.chr_id;
                   end_pos = pu.pos; 

                   long contig_len = 1000;
                    
                   pu = pr.getNext(); 
                   while (pu != null && pu.chr_id.equals(start_contig)){
                       int depth = pu.triplets.get(0).depth;
                       if (depth == min_depth)
                           break;

                       end_pos = pu.pos;
                       pu = pr.getNext();
                   }
                 
                   start_pos = tail_length;
                   
                   regions.add(new SequenceRegion(start_contig, start_pos, end_pos));
                }else
                    break; 
           }//while true
            System.out.println("Number of Regions:" + regions.size());
           
           int chain_start=0;
           while( chain_start<regions.size()){
               int i=chain_start;
               merged_regions.add(new SequenceRegion(regions.get(chain_start).contig, regions.get(chain_start).start, regions.get(i).end));
               chain_start = i+1;
           }
           System.out.println("Number of Merged Regions:" + merged_regions.size());
           
           for(SequenceRegion sr : merged_regions)
               writer.write(String.format("%s\t%d\t%d\n", sr.contig, sr.start-1, sr.end));
           
       
        
    }

    
}
