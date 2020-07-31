package tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;

public class EnvMaker {
    private  int pLens[], count;
    private  String pNames[], primers[];
    private DecimalFormat df;
    private File inFile;
    private char [] allow;
    private int seqCount;
    private int primerNum;

    public EnvMaker(String pFilename, boolean hasProbe) {
        inFile = new File(pFilename);
        this.count = 0;

        df = new DecimalFormat();
        df.setMaximumFractionDigits(2);

        seqCount=0;

        if(hasProbe){
            primerNum = 6;
        } else {
            primerNum = 4;
        }

        buildAllowList();

        primers=new String[primerNum];
        for(int i=0;i<primers.length;i++) {
            primers[i]="";
        }

        pNames=new String[primers.length];

        if(hasProbe) {
            pNames[0] = "5'-Fwd-3'-";//matches onto 3'-5' strand
            pNames[1] = "5'-Probe-3'-";//matches onto 3'-5' strand
            pNames[2] = "5'-Rev-3'-";//matches onto 3'-5' strand
            pNames[3] = "3'-Fwd-5'-";//matches onto 5'-3' strand
            pNames[4] = "3'-Probe-5'-";//matches onto 5'-3' strand
            pNames[5] = "3'-Rev-5'-";//matches onto 5'-3' strand
        } else {
            pNames[0] = "5'-Fwd-3'-";//matches onto 3'-5' strand
            pNames[1] = "5'-Rev-3'-";//matches onto 3'-5' strand
            pNames[2] = "3'-Fwd-5'-";//matches onto 5'-3' strand
            pNames[3] = "3'-Rev-5'-";//matches onto 5'-3' strand
        }

        checkNumPrimers();
    }

    public void checkNumPrimers(){
        count=0;
        try {
            BufferedReader input =  new BufferedReader(new FileReader(inFile));

            try {
                String line = null;

                while (( line = input.readLine()) != null) {
                    if(line.indexOf(">")==0) {
                        count++;

                        //Primer names
                        pNames[count-1]+=line.substring(1, line.length());
                        pNames[count-1+(primerNum / 2)]+=line.substring(1, line.length());

                        if(count>3) {
                            System.out.println("Error - more than 3 sequences in primer file");
                            System.out.println("Exiting...");
                            System.exit(0);
                        }
                    }
                    else {
                        primers[count-1]+=line;
                    }
                }
            }
            finally {
                input.close();
            }
        }
        catch (IOException ex) {
            ex.printStackTrace();
        }

        if(count==2) {
//            System.out.println("Only two sequences in primers file - creating blank probe sequence");
//            primers[2]=primers[1];
//            primers[5]=primers[1];
//            primers[1]="NNNNNNNNNNNNNNNN";
//            primers[4]=primers[1];
            System.out.println("Forward Primer Seq = "+primers[0]+" length = "+primers[0].length());
            System.out.println("Reverse Primer Seq = "+primers[1]+" length = "+primers[1].length());
        }
        else if(count!=3) {
            System.out.println("Error - abnormal number of sequences in primer file (expected 2 primers and/or 1 probe");
            System.out.println("Exiting...");
            System.exit(0);
        }
        if (count == 3) {
            System.out.println("Forward Primer Seq = " + primers[0] + " length = " + primers[0].length());
            System.out.println("Probe Seq = " + primers[1] + " length = " + primers[1].length());
            System.out.println("Reverse Primer Seq = " + primers[2] + " length = " + primers[2].length());
        }
    }

    public String [] setPrimers(){
        //CheckPrimers
        for(int p=0;p<(primerNum/2);p++) {
            //Technically GAPs are in the 'allowables' and get removed but this is for future development
            if(primers[p].indexOf("-")>=0) {
                System.out.println("Error - primer sequence has a gap in it = "+primers[p]);
                System.out.println("Exiting...");
                System.exit(0);
            }

            primers[p]=checkSeqs(pNames[p],primers[p]);
        }

        //ReversePrimers - NOT reverse-complement
        //The original 5'3' primers are checked against the complement target sequence which are in 3'-5'
        //Then check the reverse primers (now 3'-5') against the original 5'-3' target seqs
        for(int p=0;p<(primerNum/2);p++) {
            primers[p+(primerNum/2)]=reverse(primers[p]);
        }

        pLens=new int[primers.length];
        for(int p=0;p<primers.length;p++) {
            pLens[p]=primers[p].length();
        }

        return primers;
    }

    public  String checkSeqs(String name, String inSeq) {

        String outSeq="";
        String bad="";
        int b=0;

        //Switch Us to Ts, Remove GAPs
        outSeq=inSeq.toUpperCase().replace("U", "T").replace("-", "");

        for(int i=0;i<outSeq.length();i++) {
            boolean test=false;

            for(int j=0;j<allow.length;j++) {
                if(outSeq.charAt(i)==allow[j]) {
                    test=true;
                    break;
                }
            }
            if(!test) {
                bad+=outSeq.charAt(i)+" ";
                b++;
            }
        }

        if(b>0) {
            System.out.println(name+" had "+b+" bad seq bases = "+bad+" original seq = "+inSeq);
            System.out.println("Exiting...");
            System.exit(0);
        }

        return outSeq;
    }

    private String reverse(String inSeq) {

        String rev="";

        //Reverse
        for(int i=inSeq.length()-1;i>=0;i--) {
            rev=rev+inSeq.charAt(i);
        }

        return rev;
    }

    private void buildAllowList() {
        //Allowed DNA characters
        allow=new char[16];
        allow[0]='A';
        allow[1]='C';
        allow[2]='G';
        allow[3]='T';//1st thing we do is flip Us to Ts - so Us not allowed technically - so not in list
        allow[4]='M';
        allow[5]='R';
        allow[6]='W';
        allow[7]='S';
        allow[8]='Y';
        allow[9]='K';
        allow[10]='V';
        allow[11]='H';
        allow[12]='D';
        allow[13]='B';
        allow[14]='N';
        allow[15]='-';//left in for future alignment stuff - technically remove gaps at start currently
    }

    public String [] getpNames(){
        return pNames;
    }

    public void checkSeqs(String sFilename){
        File inFile = new File(sFilename);
        seqCount = 0;
        try {
            BufferedReader input =  new BufferedReader(new FileReader(inFile));

            try {
                String line = null;

                while (( line = input.readLine()) != null) {
                    if(line.indexOf(">")==0) {
                        seqCount++;
                    }
                }
            }
            finally {
                input.close();
            }
        }
        catch (IOException ex) {
            ex.printStackTrace();
        }

        if(seqCount==0) {
            System.out.println("Error - no seqs found - check file is FASTA format?");
            System.out.println("Exiting...");
            System.exit(0);
        }
    }
    public int getSeqCount(){
        return seqCount;
    }

    public int[] getpLens() {
        return pLens;
    }
}
