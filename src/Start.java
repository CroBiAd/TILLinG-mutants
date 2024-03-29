
import htsjdk.common.BAMFilter;
import htsjdk.common.BAMRepairer;
import htsjdk.common.MUTProcessor;
import java.io.File;
import java.io.IOException;
import java.io.StringReader;
import java.util.Arrays;
import java.util.Properties;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author Mario Fruzangohar
 */
public class Start extends javax.swing.JFrame {

    /**
     * Creates new form Start
     */
    public Start() {
        initComponents();
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 400, Short.MAX_VALUE)
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 300, Short.MAX_VALUE)
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {
        /* Set the Nimbus look and feel */
        //<editor-fold defaultstate="collapsed" desc=" Look and feel setting code (optional) ">
        /* If Nimbus (introduced in Java SE 6) is not available, stay with the default look and feel.
         * For details see http://download.oracle.com/javase/tutorial/uiswing/lookandfeel/plaf.html 
         */
        try {
            for (javax.swing.UIManager.LookAndFeelInfo info : javax.swing.UIManager.getInstalledLookAndFeels()) {
                if ("Nimbus".equals(info.getName())) {
                    javax.swing.UIManager.setLookAndFeel(info.getClassName());
                    break;
                }
            }
        } catch (ClassNotFoundException ex) {
            java.util.logging.Logger.getLogger(Start.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (InstantiationException ex) {
            java.util.logging.Logger.getLogger(Start.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (IllegalAccessException ex) {
            java.util.logging.Logger.getLogger(Start.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (javax.swing.UnsupportedLookAndFeelException ex) {
            java.util.logging.Logger.getLogger(Start.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }
        //</editor-fold>
        boolean showForm = true;
        String srcDir = null;
        String destDir = null;
        String task = null;
        String host = null; // db host
        String exporterClass =null;
        String p1 =null;  // this is used for a general purpose parameter
        String p2 =null;  // this is used for a general purpose parameter
        String p3 =null;  // this is used for a general purpose parameter
        String p4 =null;  // this is used for a general purpose parameter
        String p5 =null;  // this is used for a general purpose parameter
        String p6 =null;  // this is used for a general purpose parameter
        String p7 =null;  // this is used for a general purpose parameter
        String p8 =null;  // this is used for a general purpose parameter
        String threads =null;  // this is used for number of threads in thread pool
        
        
        String mate1File = null;
        String mate2File = null;
        String unpairedFile = null;
        
        String file1 = null;
        String file2 = null;
        String file3 = null;
        String file4 = null;
        String file5 = null;
        
        
        
        Properties prop = new Properties();
        StringBuffer sb = new StringBuffer();
        for(int i =0; i<args.length; i++){
            sb.append(args[i]+"\n");
        }
        try {
            prop.load(new StringReader(sb.toString()));
        } catch (IOException ex) {
             System.out.println(ex.getMessage());
        }

       //--------------------------------------------
        try{
           if(prop.getProperty("showform").equals("no")){
             showForm = false;
           }
           //-------------------------------------------
           task = prop.getProperty("task");
           host = prop.getProperty("host");
           srcDir = prop.getProperty("srcdir");
           destDir = prop.getProperty("destdir");
           exporterClass = prop.getProperty("class");
           p1 = prop.getProperty("p1");
           p2 = prop.getProperty("p2");
           p3 = prop.getProperty("p3");
           p4 = prop.getProperty("p4");
           p5 = prop.getProperty("p5");
           p6 = prop.getProperty("p6");
           p7 = prop.getProperty("p7");
           p8 = prop.getProperty("p8");
           threads = prop.getProperty("threads");
           
           mate1File = prop.getProperty("1");
           mate2File = prop.getProperty("2");
           unpairedFile = prop.getProperty("U");
           file1 = prop.getProperty("file1");
           file2 = prop.getProperty("file2");
           file3 = prop.getProperty("file3");
           file4 = prop.getProperty("file4");
           file5 = prop.getProperty("file5");
           

        }catch(NullPointerException ex){}
     if(showForm){
        /* Create and display the form */
        java.awt.EventQueue.invokeLater(new Runnable() {
            
            public void run() {
                new Start().setVisible(true);
            }
        });
     }else{
         if (task.equals("mutation2")){
            if(srcDir == null || destDir == null || p1 == null || p2 == null || p3 == null ){
                   System.out.println("'All following parameters are required: 'srcdir', 'destdir', 'p1', 'p2', 'p3'");  
                   System.exit(1);                     
                }
                 //developed 31/5/2018 for new sagov tilling population mutation calling
                 //This task creates 2 files unique_mutations_pr.mut and non_unique_mutations_pr.mut, the first file has following 12 columns:
                 //chr_id, pos_start, pos_end, sample_id, mutation_type, mutation_subtype,  affected_base, offset, prob, zygosity, read_depth, ref_base
                 //parameters
                 //srcDir:pileup file
                 //destdir: directory of generated mut files
                 //p1: number of samples including parent
                 //p2: min number of reads for a non-reference allele
                 //p3 : Minor allele ratio to determin zygosity can be like 0.2 ...  
                 
                
                PileupProcessor pp = new PileupProcessor();
                File dataDir = new File(srcDir);
                if (dataDir.isFile()){
                    System.out.println("Processing File: " + dataDir.getName());
                    try{
                       pp.callMutation2(dataDir.getAbsolutePath(), destDir, Integer.parseInt(p1), Integer.parseInt(p2), Double.parseDouble(p3), 1);
                    }catch(Exception ex){}
                }
                System.out.println("Done.");

     }else if (task.equals("buildsampleseq")){
                //Developed on 5/6/2018
                // This task accepts one mutation file and produces a new mutation file
                // This task creates a new mutation file with 1 extra columns, this extra column is reference replaced by sample genotype at each position, 
                //it does not try resolve haplotype in polymorphic positions, it just report polymorphic position in square bracket such as [A-] [+AG] mutation itself is also in a bracket [A/-] 
                //mutation file accepted by this task is 8 column format 
                
                if(srcDir == null  || p1==null || p2 ==null  || mate1File==null || mate2File==null || file1==null){
                   System.out.println("'All following parameters are required: 'srcdir', '1', '2', 'p1' , 'p2' , 'file1' ");  
                   System.exit(1);                     
                }
                //src dir = mutation file                
                //p1 = min flanking mutation position
                //p2 = max flanking mutation position
              
                
                //1=fasta reference file
                //2=fasta index file
                //file1=genotype file that only shows differences with reference (the output of genotype2 task)
                System.out.println("processing mutation file  " + srcDir);
                MUTProcessor processor = new MUTProcessor(srcDir, true, true); 
                processor.buildSampleSeq(Integer.parseInt(p1), Integer.parseInt(p2), mate1File, mate2File, file1.trim()); 
                 

                System.out.println("Done.");
            }else if (task.equals("genotype2")){
                if(srcDir == null || destDir == null || p1 == null || p2 == null  ){
                   System.out.println("'All following parameters are required: 'srcdir', 'destdir', 'p1', 'p2'");  
                   System.exit(1);                     
                }
                //This is for my own Durum Diversity Panel project
                 //parameters
                 //srcDir: pileup file
                 //destdir: generated genotype file (Tab separated)
                 //p1: number of samples
                 //p2: min number of reads for a non-reference allele
                 //p3 (optional):use MAPQ if it is set, otherwise not using mapq
                 //p4 :Minor allele ratio to determin zygosity can be like 0.2 ...
                 //p5 (optional): output types : 1)All   2)differences (see the function callGenotype2 for definition)
                boolean all_positions = true;
                boolean useMAPQ = false;
                if (p3!=null)
                    useMAPQ = true;
                if (p5!=null){
                    if (Byte.parseByte(p5) == 2)
                        all_positions = false; // then we just want to output positions that are different from the reference
                }
                
                htsjdk.common.PileupProcessor pp = new htsjdk.common.PileupProcessor();
                File dataDir = new File(srcDir);
                if (dataDir.isFile()){
                    System.out.println("Processing File: " + dataDir.getName());                    
                    pp.callGenotype2(dataDir.getAbsolutePath(), destDir, Integer.parseInt(p1), Integer.parseInt(p2), Double.parseDouble(p4), useMAPQ, all_positions);  
                }
                System.out.println("Done.");

            }else if (task.equals("bamrepair")){
                if(srcDir == null || destDir == null){
                   System.out.println("'All following parameters are required: 'srcdir', 'destDir'");  
                   System.exit(1);                     
                }
                BAMRepairer bp = new BAMRepairer();
                File dataDir = new File(srcDir);
                if (dataDir.isFile()){
                    System.out.println("Processing File: " + dataDir.getName());
                    File f = new File(destDir, dataDir.getName());
                    bp.remove_duplicate_reads(dataDir.getAbsolutePath(), f.getAbsolutePath() );
                }else{                 
                    File[] files = dataDir.listFiles(new BAMFilter());
                    Arrays.sort(files);

                    for (int i = 0; i < files.length; ++i){
                        if(files[i].isFile()){
                            File f = new File(destDir, files[i].getName());
                            System.out.println("Processing File: " + f.getName());
                            bp.remove_duplicate_reads(files[i].getAbsolutePath(), f.getAbsolutePath() );
                            //System.out.println("fastq file=" + files[i].getAbsolutePath());
                        }

                    }
                }
                System.out.println("Done.");
            }else if (task.equals("highcoverage")){
                if(file1== null || srcDir == null || destDir == null || p1==null || p2==null || p3==null){
                   System.out.println("'All following parameters are required: 'file1', srcdir', 'destdir', 'p1', 'p2', 'p3'");  
                   System.exit(1);                     
                }  
                //file1: fasta index file
                //srcdir: pileup file
                //destdir : produced bed files
                //p1: minimum depth
                //p2: tail length
                //p3: maximum distance to merge 2 coninuous regions
                
                PileupProcessor pp = new PileupProcessor();
                File pileup = new File(srcDir);
                if (pileup.isFile()){
                    System.out.println("Processing File: " + pileup.getName());
                    pp.disocoverHighCoverage(pileup.getAbsolutePath(), file1 , Integer.parseInt(p1), Integer.parseInt(p2), Integer.parseInt(p3), destDir);     
                }
                System.out.println("Done.");
              
            }
    }
}
    // Variables declaration - do not modify//GEN-BEGIN:variables
    // End of variables declaration//GEN-END:variables
}
