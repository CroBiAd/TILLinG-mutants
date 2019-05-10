
import htsjdk.samtools.BAMRecord;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import java.io.File;
import java.io.IOException;
import java.util.Iterator;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author Mario Fruzangohar
 */
public class BAMRepairer {
         public void remove_duplicate_reads(String bamFile, String outFile){
       
        SamReaderFactory factory = SamReaderFactory.makeDefault().enable(SamReaderFactory.Option.DONT_MEMORY_MAP_INDEX);
        SamReader reader = factory.open(new File(bamFile));
        FileWriter writer = null; 
        int rec_count = 0; 
        int sm_count = 0; 
        int pm_count = 0;           
        int q0_count = 0;
        int q1_count = 0;
        
         try{
           Iterator it = reader.iterator();
           boolean isPairedMapped = false;
           
           boolean last_strand = false;           
           String last_seq = "";
           String last_contig = "";
           int last_pos = 0;
           
           String last_mate_contig = "";
           int last_mate_pos = 0;
           boolean last_mate_strand = false;
           // because bam file is sorted, then alignment are sorted based on their position on the contig 
           while (it.hasNext()){
               rec_count++;
               BAMRecord rec = it.next();
              
               boolean strand = rec.getReadNegativeStrandFlag();
               
               String contig = rec.getReferenceName();
               int pos = rec.getAlignmentStart(); 
               String seq = rec.getReadString();

               isPairedMapped = false;
               String mate_contig = "";
           

               if (rec.getReadPairedFlag() ){
                   isPairedMapped = true;
                                      
                   mate_pos = rec.getMateAlignmentStart();
                   mate_strand = rec.getMateNegativeStrandFlag();                   
               }
              
               boolean A = contig.equals(last_contig) ;
               boolean B = mate_contig.equals(last_mate_contig) ;
               
               if (!A || !B){
                       writer.write(rec);                   
                       last_strand = strand;
                       
                       last_seq = seq;
                       last_contig = contig;
                       last_pos = pos;
                       
                       
                       last_mate_contig = mate_contig;
                 
                       
                   
               }else{
                   
                   if (rec.getMappingQuality() == 0)
                           q0_count++;
                   if (rec.getMappingQuality() == 1)
                           q1_count++;

                   if (!isPairedMapped)
                       sm_count++; 
                   else
                       pm_count++; 
                   
               }// else
           }//while records
           System.out.println(String.format("Total Number of Records: %,d", rec_count) );
           System.out.println(String.format("Total Number of Single Mapped (not single-end) Reads Removed, : %,d  Paired Mapped Reads Removed: %,d  (SUM:%,d) Quality Zero removed: %,d   Quality One removed: %,d", sm_count, pm_count, sm_count+pm_count, q0_count, q1_count) );           
       }finally{
           try {
               reader.close();
               writer.close();
           } catch (IOException ex) {}
       }
    }

}
