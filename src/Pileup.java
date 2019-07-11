
import htsjdk.common.Mutation;
import htsjdk.common.PileupTriplet;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author Mario Fruzangohar
 */
public class Pileup extends htsjdk.common.Pileup{
        public List<Mutation> callMutationsPerSample2(int mutant_sample_count, int non_ref_count, double minor_ratio) throws Exception{
        List<Mutation> mutations = new ArrayList();
        PileupTriplet ctrl = triplets.get(0);
        ctrl.countAll(false);      
        for(int i=1; i<=mutant_sample_count; i++)
            triplets.get(i).countAll(false);      
          
        Set<Byte> ctrlBaseSet = ctrl.getGenotypeSet(non_ref_count); 
              
        if(ctrl.depth == 0) 
            return mutations; 
  
        
       
        Map<String, Integer> ctrlInsertMap = ctrl.insertmap;        

        
       
        for(int i=1; i<=mutant_sample_count; i++){
            PileupTriplet sa = triplets.get(i);
            String mutant_gstr = sa.getDiploidGenotype( minor_ratio, non_ref_count, false);
            
            if (sa.depth == 0)
                continue; 
            Set<Byte> mBaseSet  = sa.getGenotypeSet(non_ref_count);  
            mBaseSet.removeAll(ctrlBaseSet);
            if(mBaseSet.size()>0){
                for(byte mutation_allele: mBaseSet){
                    for(byte ctrl_allele : ctrlBaseSet)
                        if(ctrl_allele!=mutation_allele && ctrl_allele!=htsjdk.common.Helper.D){
                           Mutation mu = new Mutation();
                           mu.sample_id = i ;
                           mu.refBase = sa.parent.ref_nt;
                           mu.oldBase = ctrl_allele;
                           if (mutant_gstr.length()>1)
                               mu.zygosity = htsjdk.common.Helper.ZG_HETROZYGOUS;
                           else
                               mu.zygosity = htsjdk.common.Helper.ZG_HOMOZYGOUS;
                           
                           if(mutation_allele==htsjdk.common.Helper.D){
                               mu.type = htsjdk.common.Helper.MT_DELETION; 
                               mu.read_depth = sa.deleted;
                           }else{
                               mu.type = htsjdk.common.Helper.MT_Substitution;
                               mu.newBase = mutation_allele;
                               mu.read_depth = sa.counts[mutation_allele];
                           }
                           mutations.add(mu);

                        }
                }                
            }
            
        }  // for each sample
        
        
        return mutations;
    }

}
