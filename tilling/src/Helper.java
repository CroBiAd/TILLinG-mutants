
import static htsjdk.common.Helper.A;
import static htsjdk.common.Helper.C;
import static htsjdk.common.Helper.D;
import static htsjdk.common.Helper.G;
import static htsjdk.common.Helper.GAP;
import static htsjdk.common.Helper.MT_DELETION;
import static htsjdk.common.Helper.MT_Default;
import static htsjdk.common.Helper.MT_Hidden;
import static htsjdk.common.Helper.MT_INSERTION;
import static htsjdk.common.Helper.MT_Substitution;
import static htsjdk.common.Helper.N;
import static htsjdk.common.Helper.T;
import static htsjdk.common.Helper.UNKNOWN;
import static htsjdk.common.Helper.ZG_HETROZYGOUS;
import static htsjdk.common.Helper.ZG_HOMOZYGOUS;
import static htsjdk.common.Helper.ZG_INVALID;
import static htsjdk.common.Helper.ZG_UNKNOWN;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author Mario Fruzangohar
 */
public class Helper {
      // number neucleotie letters alphabetically
   public static final byte A = 0;   
   public static final byte C = 1;
   public static final byte G = 2;
   public static final byte T = 3;
   public static final byte GAP = 5; // in pairwise and mutiple alignment this represent a GAP in the alignment
   public static final byte N = 6; // any base
   public static final byte UNKNOWN = -1;
   public static final byte D = -2; // base is deleted
   
  public static final byte ZG_UNKNOWN = 0; // Zygosity unknown
   public static final byte ZG_HOMOZYGOUS = 1; // Homozygous
   public static final byte ZG_HETROZYGOUS = 2; // Hetrozygous   
   public static final byte ZG_INVALID = -1; // something is not right!   
   
 public static final byte MT_DELETION = 1; // base  deletion
   public static final byte MT_INSERTION = 2; // base  insertion
   public static final byte MT_Substitution = 3; // base  change
   public static final byte MT_Default = 4; // default mutation type
   public static final byte MT_Hidden = 5; // This is subtype of Transition mutation, when ratio test will reveal mutation
    public static String Base_String(byte b){
       String ch = " ";
       switch(b){
           case A:
               ch = "A";
               break;
           case T:
               ch = "T";
               break;
           case C:
               ch = "C";
               break;
           case G:
               ch = "G";
               break;
            case N:
               ch = "N";
               break; 
            case GAP:
                ch = "-";
                break;
            case D:
                ch = "D";
                break;

            case UNKNOWN:
                ch = "";
                break;
                
       }
       
       return ch;
   }
    public static String MT_Short_String(byte mt){
       String str = "";
       switch(mt){
           case MT_DELETION:
               str = "DEL";
               break;
           case MT_INSERTION:
               str = "INS";
               break;
           case MT_Substitution:
               str = "SUB";
               break;
           case MT_Default:
               str = "DEF";
               break;
           case MT_Hidden:
               str = "HID";
               break;
               
       }
       return str;
   }
    
     public static byte nt_to_i(char nt){
           switch (nt){
                case 'A':
                case 'a':
                    return A;                       
                case 'T':                        
                case 't':
                    return T;                        
                case 'C':
                case 'c':
                    return C;                        
                case 'G':                        
                case 'g':
                    return G;   
                case '-':
                    return GAP;
                case 'n':
                case 'N':
                    return N;
                case 'd':
                case 'D':
                case '*':
                    return D; 
            }
           return UNKNOWN;       
    }  
     public static String ZG_Short_String(byte zg){
       String str = "";
       switch(zg){
           case ZG_UNKNOWN:
               str = "UNKNOWN";
               break;
           case ZG_HOMOZYGOUS:
               str = "HOM";
               break;
           case ZG_HETROZYGOUS:
               str = "HET";
               break;
           case ZG_INVALID:
               str = "ERR";
               break;           
       }
       return str;
   }
     
     public static boolean isVariationFromReference(String[] genotypes, String ref){
              for (String g :genotypes)
           if(g!=null && !g.isEmpty() &&  !g.toLowerCase().equals(ref.toLowerCase()))
              return true; 
       
       
        return false;  
   }
     
     
}
