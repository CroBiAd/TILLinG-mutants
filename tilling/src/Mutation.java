/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author Mario Fruzangohar
 */
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */


import java.util.ArrayList;
import java.util.List;
import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.commons.lang.builder.HashCodeBuilder;

/**
 *
 * @author mario.fruzangohar
 * 
 */
public class Mutation extends htsjdk.common.Mutation{
 
    
    public boolean unique = true; 
    public int sample_id = -1; 
    public String sample_name = null; 
    public int offset = 0; 
    public byte type = Helper.MT_Default;
    public byte subtype = Helper.MT_Default; 
    public byte zygosity = Helper.ZG_UNKNOWN; 
    public byte refBase = -1;
    public byte oldBase = -1;
    public byte newBase = -1;
    public String insertStr = "";
    
    
    public double prob = 0; 
    public int read_depth = 0; 
   
    public String contig = null;
    public int position = 0; 
    
    public boolean equals(Object other){
        if (other == null) return false;
        if (other == this) return true;
        if (!(other instanceof Mutation))return false;
        Mutation mu = (Mutation)other;
        return new EqualsBuilder().append(type, mu.type).append(subtype, mu.subtype).append(oldBase, mu.oldBase).append(newBase, mu.newBase).append(insertStr, mu.insertStr).append(offset, mu.offset).isEquals();
        
    }
    
    
    public String print12(){
        //12 columns format
        String out = "";
        String affected_base = "";
        String ref_base = "";
        if (type == Helper.MT_Substitution)
            affected_base = Helper.Base_String(oldBase)+Helper.Base_String(newBase);
        else if (type == Helper.MT_INSERTION)
            affected_base = insertStr ;
        else if (type == Helper.MT_DELETION)
            affected_base = Helper.Base_String(oldBase);
        
        if(refBase!=-1)
            ref_base = Helper.Base_String(refBase);
        if (unique){            
              out = String.format("%s\t%s\t%s\t%s\t%d\t%.4f\t%s\t%s", sample_name, Helper.MT_Short_String(type), Helper.MT_Short_String(subtype),  affected_base, offset, prob, Helper.ZG_Short_String(zygosity), ref_base);
        }else{
            
           out = String.format("%s\t%s\t%s\t%d\t%.4f\t%s", Helper.MT_Short_String(type), Helper.MT_Short_String(subtype), affected_base, offset, prob, ref_base);
        }
        
        return out;
    }
    
    
    public void getFromString8(String str){
        //contig20766_1_2204_3A_10388	2046	G	453	SUB	GA	HET	355
        String[] arr = str.split("\\s");
        contig = arr[0].trim();
        position = Integer.parseInt(arr[1]);
        refBase = Helper.nt_to_i(arr[2].charAt(0));
        sample_id = Integer.parseInt(arr[3]);
        if (arr[4].startsWith("TRN") || arr[4].startsWith("SUB"))
            type = Helper.MT_Substitution;
        else if (arr[4].startsWith("INS"))
            type = Helper.MT_INSERTION;
        else if (arr[4].startsWith("DEL"))
            type = Helper.MT_DELETION;

        if (type == Helper.MT_Substitution){
            oldBase = Helper.nt_to_i(arr[5].charAt(0)) ;
            newBase = Helper.nt_to_i(arr[5].charAt(1)) ;
        }else if (type == Helper.MT_INSERTION){
            insertStr = arr[5];
        }else if (type == Helper.MT_DELETION){
            oldBase = Helper.nt_to_i(arr[5].charAt(0)) ;
        }
       if(arr[6].startsWith("HOM"))
            zygosity = Helper.ZG_HOMOZYGOUS;
        else if (arr[6].startsWith("HET"))
            zygosity = Helper.ZG_HETROZYGOUS;
        else if (arr[6].startsWith("UNKNOWN"))
            zygosity = Helper.ZG_UNKNOWN;        
        //sample_name = arr[3];
       
       read_depth = Integer.parseInt(arr[7]);
    }
    public void getFromString12(String str){
        // this works opposite of print 12 columns format
        String[] arr = str.split("\\s");
        contig = arr[0].trim();
        position = Integer.parseInt(arr[1]);
        //sample_id = Integer.parseInt(arr[3]);
        sample_name = arr[3];
        if (arr[4].startsWith("TRN") || arr[4].startsWith("SUB"))
            type = Helper.MT_Substitution;
        else if (arr[4].startsWith("INS"))
            type = Helper.MT_INSERTION;
        else if (arr[4].startsWith("DEL"))
            type = Helper.MT_DELETION;
        
        if (arr[5].startsWith("DEF"))
            subtype = Helper.MT_Default;
        else if (arr[5].startsWith("INS"))
            subtype = Helper.MT_INSERTION;
        else if (arr[5].startsWith("DEL"))
            subtype = Helper.MT_DELETION;
        else if (arr[5].startsWith("HID"))
            subtype = Helper.MT_Hidden;
        
        if (type == Helper.MT_Substitution){
            oldBase = Helper.nt_to_i(arr[6].charAt(0)) ;
            newBase = Helper.nt_to_i(arr[6].charAt(1)) ;
        }else if (type == Helper.MT_INSERTION){
            insertStr = arr[6];
        }else if (type == Helper.MT_DELETION){
            oldBase = Helper.nt_to_i(arr[6].charAt(0)) ;
        }

        offset = Integer.parseInt(arr[7]);
        prob = Double.parseDouble(arr[8]);
        
        if(arr[9].startsWith("HOM"))
            zygosity = Helper.ZG_HOMOZYGOUS;
        else if (arr[9].startsWith("HET"))
            zygosity = Helper.ZG_HETROZYGOUS;
        else if (arr[9].startsWith("UNKNOWN"))
            zygosity = Helper.ZG_UNKNOWN;
        
        if (arr[10].length()>0)
            refBase = Helper.nt_to_i(arr[10].charAt(0)); 
        
        read_depth = Integer.parseInt(arr[11]);
    }
    
    public String print2(){
        //2 columns format
        String out = "";
        String affected_base = "";       
        if (type == Helper.MT_Substitution)
            affected_base = Helper.Base_String(oldBase)+Helper.Base_String(newBase);
        else if (type == Helper.MT_INSERTION)
            affected_base = insertStr ;
        else if (type == Helper.MT_DELETION)
            affected_base = Helper.Base_String(oldBase);
        
        
            
        out = String.format("%s\t%s", Helper.MT_Short_String(type), affected_base);
        
        
        return out;
    }
    
    public String print4website(){
        String out = null;
       if (type == Helper.MT_Substitution)
           out = String.format("[%s/%s]", Helper.Base_String(oldBase), Helper.Base_String(newBase));
        else if (type == Helper.MT_INSERTION)
            out = String.format("[-/%s]",  insertStr);
        else if (type == Helper.MT_DELETION)
            out = String.format("[%s/-]", Helper.Base_String(oldBase));
        
        
        return out; 
    }
    
}
