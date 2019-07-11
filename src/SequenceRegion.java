/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author Mario Fruzangohar
 */
public class SequenceRegion {
   public String contig = null;
    public long start = 0;  //start is always less than end
    public long end = 0;
    
    public SequenceRegion(){
        
    }
    
    public SequenceRegion(String cntg, long s, long e){
        this.contig = cntg;
        this.start = s;
        this.end = e;        
    }   
}
