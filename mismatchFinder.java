//EDIT of GoPrime.java to run with or without Probe and give only mismatch info

import tools.EnvMaker;

import java.io.*;
import java.text.DecimalFormat;

public class misMatchFinder {

    public int primerNum;
    public  String primers[],  seqs[], compSeqs[], seqNames[], pNames[];
    public  String outFilename, head, header, tTx;
    public  static String pFilename, sFilename;

    public  int scores[][][], pLens[], mLen, count;

    public  double minPm, minPmI, minPb;

    public  char bases[][], revCodes[][],  amb[][], comps[][], compsBase[], ts[][], muts[][];

    public  double max14nt;

    public  static DecimalFormat df;

    public  boolean success, hasProbe;

    public  String noHits=""+System.getProperty("line.separator");



    public static void main(String[] args) {

        System.out.println("GoPrime Started" + System.getProperty("line.separator"));

        df = new DecimalFormat();
        df.setMaximumFractionDigits(2);
        boolean has_Probe = false;
        if (args.length == 2) {
            pFilename = args[0];
            sFilename = args[1];
            has_Probe = false;
        } else if (args.length ==3){
            pFilename = args[0];
            sFilename = args[1];
            has_Probe = true;
        } else {
            System.out.println("Error - incorrect usage: " + args.length + " parameters were given");
            System.out.println("Correct usage = java -jar GoPrime.jar primerSeqFilename.fasta targetSeqFilename.fasta");
            System.out.println("Exiting..." + System.getProperty("line.separator"));
            System.exit(0);
        }

        System.out.println("Primer Sequence File = " + pFilename);
        System.out.println("Target Sequence File = " + sFilename);

        //Initiate
        misMatchFinder finder = new misMatchFinder();
        finder.Run(has_Probe);
    }
    public void Run(boolean has_probe){
        //PRINT Output file names
        hasProbe = has_probe;
        if (hasProbe){
            primerNum = 6;
        } else {
            primerNum = 4;
        }

        setUpBases();
        setUpCts();

        EnvMaker envMaker = new EnvMaker(pFilename, hasProbe);

        primers = envMaker.setPrimers();
        pNames = envMaker.getpNames();
        pLens = envMaker.getpLens();

        //Combined primer length - needed for combined mismatch freq
        mLen = primers[0].length() + primers[primerNum /2 -1].length();

        envMaker.checkSeqs(sFilename);

        buildSeqList(envMaker, envMaker.getSeqCount());
        //End Setup

        try {
            outFilename = sFilename + "_out.txt";

            FileWriter fstreamOut = new FileWriter(outFilename);
            BufferedWriter outOut = new BufferedWriter(fstreamOut);

            outFilename = sFilename + "_cts.txt";

            FileWriter fstreamCT = new FileWriter(outFilename);
            BufferedWriter outCT = new BufferedWriter(fstreamCT);

            //outputPrinter printer = new outputPrinter(sFilename, pLens);

            //Loop through each sequence in turn
            for (int s = 0; s < compSeqs.length; s++) {
                success = false;

                try {

                    //outFilename=sFilename.substring(0, sFilename.lastIndexOf("."))+"_"+s+"_matrix.txt";
                    outFilename=sFilename+"_"+(s+1)+"_matrix.txt";

                    tTx="Evaluating seq "+(s+1)+" "+seqNames[s];
                    System.out.println(tTx);
                    outOut.write(tTx+System.getProperty("line.separator"));

                    //Clear results - 6 for F/P/R*2, 14 for all the different mutation counts - plus one for Ns now
                    scores = new int[compSeqs[s].length()][primerNum][16];

                    runMismatchDetection(s);

                    PrintCandidatePos(outOut);
                }//matrix outputFile try

                catch (Exception e) {
                    e.printStackTrace();
                    System.err.println("Error: " + e.getMessage());
                }

                if(!success) {
                    noHits+=seqNames[s]+"\tNOHIT"+System.getProperty("line.separator");
                }
            }//end of comSeq loop s - evaluate each seq

            outOut.close();

            outCT.write(noHits);
            outCT.close();
        }//cts outFile try
        catch (Exception e) {
            e.printStackTrace();
            System.err.println("Error: " + e.getMessage());
        }
        System.out.println("GoPrime finished...Exiting");
    }

    private  int checkMutation(char pBase, char tBase) {

        int mut=0;
        int pSel=0, tSel=0;

        for(int i=0;i<compsBase.length;i++) {
            if(pBase==compsBase[i]) {
                pSel=i;
                break;
            }
        }

        for(int i=0;i<compsBase.length;i++) {
            if(tBase==compsBase[i]) {
                tSel=i;
                break;
            }
        }

        //check if it is a match
        boolean test=false;
        for(int i=0;i<comps[pSel].length;i++) {
            if(comps[pSel][i]==tBase) {
                test=true;
                break;
            }
        }

        //if not - check if the mutation is a transition or transversion
        if(!test) {
            //New - we have to revcomp the primer base
            //This is what we expect in the template
            //Compare this to what we observe to determine if it is mutated or not
            pBase=compBase(pBase);

            for(int i=0;i<compsBase.length;i++) {
                if(pBase==compsBase[i]) {
                    pSel=i;
                    break;
                }
            }

            boolean transi=false;
            for(int i=1;i<muts[pSel].length;i++) {
                if(transi)
                    break;

                for(int j=1;j<muts[tSel].length;j++) {
                    if((muts[pSel][i]==ts[0][0] & muts[tSel][j]==ts[0][1])
                            |(muts[pSel][i]==ts[0][1] & muts[tSel][j]==ts[0][0])
                            |(muts[pSel][i]==ts[1][0] & muts[tSel][j]==ts[1][1])
                            |(muts[pSel][i]==ts[1][1] & muts[tSel][j]==ts[1][0])) {
                        mut=1;
                        transi=true;
                        break;
                    }
                }
            }
            if(!transi) {
                mut=2;
            }

        }

        return mut;
    }

    private String reverseComplement(String inSeq) {

        String rev=reverse(inSeq);
        String comp=complement(rev);

        return comp;
    }

    private String reverse(String inSeq) {

        String rev="";

        //Reverse
        for(int i=inSeq.length()-1;i>=0;i--) {
            rev=rev+inSeq.charAt(i);
        }

        return rev;
    }


    private char compBase(char inBase) {

        //Complement
        char comp='N';

        if(inBase=='A') {
            comp='T';
        }
        else if(inBase=='C') {
            comp='G';
        }
        else if(inBase=='G') {
            comp='C';
        }
        else if(inBase=='T') {
            comp='A';
        }
        else {
            boolean test=false;
            for(int j=0;j<revCodes.length;j++) {
                if(revCodes[j][0]==inBase) {
                    comp=revCodes[j][1];
                    test=true;
                    break;
                }
            }
            if(!test) {
                System.out.println("ERROR - unrecognised character during complementation "+inBase);
            }
        }

        return comp;
    }


    private String complement(String inSeq) {

        //Complement
        String comp="";
        for(int i=0;i<inSeq.length();i++) {
            if(inSeq.charAt(i)=='A') {
                comp=comp+"T";
            }
            else if(inSeq.charAt(i)=='C') {
                comp=comp+"G";
            }
            else if(inSeq.charAt(i)=='G') {
                comp=comp+"C";
            }
            else if(inSeq.charAt(i)=='T') {
                comp=comp+"A";
            }
            else {
                boolean test=false;
                for(int j=0;j<revCodes.length;j++) {
                    if(revCodes[j][0]==inSeq.charAt(i)) {
                        comp=comp+revCodes[j][1];
                        test=true;
                        break;
                    }
                }
                if(!test) {
                    System.out.println("ERROR - unrecognised character during complementation "+inSeq.charAt(i));
                }
            }
        }

        return comp;
    }



    private  void setUpBases() {

        head="TotMis\tTotTs\tTotTv\tNt1Ts\tNt1Tv\tNt2Ts\tNt2Tv\tNt3Ts\tNt3Tv\tNt4Ts\tNt4Tv\t1-4NtMis\tMisPass\tPass\tNs\t5-1-4NtMis";
        header="";
        header="Pos\t"+"[Fwd]"+head+"\t[Probe]"+head+"\t[Rev]"+head+"\t[Fwd]"+head+"\t[Probe]"+head+"\t[Rev]"+head+System.getProperty("line.separator");


        //Bases and revComps
        bases=new char[4][2];
        bases[0][0]='A';bases[0][1]='T';
        bases[1][0]='C';bases[1][1]='G';
        bases[2][0]='G';bases[2][1]='C';
        bases[3][0]='T';bases[3][1]='A';

        //The Reverse Complement of IUPAC Ambiguity codes
        revCodes=new char[12][2];
        revCodes[0][0]='M';revCodes[0][1]='K';
        revCodes[1][0]='R';revCodes[1][1]='Y';
        revCodes[2][0]='W';revCodes[2][1]='W';
        revCodes[3][0]='S';revCodes[3][1]='S';
        revCodes[4][0]='Y';revCodes[4][1]='R';
        revCodes[5][0]='K';revCodes[5][1]='M';
        revCodes[6][0]='V';revCodes[6][1]='B';
        revCodes[7][0]='H';revCodes[7][1]='D';
        revCodes[8][0]='D';revCodes[8][1]='H';
        revCodes[9][0]='B';revCodes[9][1]='V';
        revCodes[10][0]='N';revCodes[10][1]='N';
        revCodes[11][0]='-';revCodes[11][1]='-';

        //This is what IUPAC amigiuity codes correspond to
        amb=new char[11][];
        amb[0]=new char[3];
        amb[0][0]='M';amb[0][1]='A';amb[0][2]='C';
        amb[1]=new char[3];
        amb[1][0]='R';amb[1][1]='A';amb[1][2]='G';
        amb[2]=new char[3];
        amb[2][0]='W';amb[2][1]='A';amb[2][2]='T';
        amb[3]=new char[3];
        amb[3][0]='S';amb[3][1]='C';amb[3][2]='G';
        amb[4]=new char[3];
        amb[4][0]='Y';amb[4][1]='C';amb[4][2]='T';
        amb[5]=new char[3];
        amb[5][0]='K';amb[5][1]='G';amb[5][2]='T';
        amb[6]=new char[4];
        amb[6][0]='V';amb[6][1]='A';amb[6][2]='C';amb[6][3]='G';
        amb[7]=new char[4];
        amb[7][0]='H';amb[7][1]='A';amb[7][2]='C';amb[7][3]='T';
        amb[8]=new char[4];
        amb[8][0]='D';amb[8][1]='A';amb[8][2]='G';amb[8][3]='T';
        amb[9]=new char[4];
        amb[9][0]='B';amb[9][1]='C';amb[9][2]='G';amb[9][3]='T';
        amb[10]=new char[5];
        amb[10][0]='N';amb[10][1]='A';amb[10][2]='C';amb[10][3]='G';amb[10][4]='T';

        int cCount=0;
        String cText="";
        char bComps[][]=new char[4][];

        //loop to get all the reverse complements of each base including ambiguity codes
        for(int b=0;b<bases.length;b++) {
            cCount=1;
            cText=""+bases[b][1];

            for(int i=0;i<amb.length;i++) {
                for(int j=1;j<amb[i].length;j++) {
                    if(amb[i][j]==bases[b][1]) {//NB - [b][1] - searching for the reverse complement of the base [b][0]
                        cText+=amb[i][0];
                        cCount++;
                    }
                }
            }

            //for each base b - store all possible revcomps bases/amiguities
            bComps[b]=new char[cText.length()];
            for(int i=0;i<cText.length();i++) {
                bComps[b][i]=cText.charAt(i);
            }
        }

        char aComps[][]=new char[amb.length][];

        //loop to find all the reverse complements of each amiguity code
        for(int a=0;a<amb.length;a++) {
            cCount=0;
            cText="";

            //loop through all the bases the amiguity stands for
            for(int j=1;j<amb[a].length;j++) {
                //for each base that the amiguity represents - find all its complements (incl ambiguity codes)
                for(int k=0;k<bases.length;k++) {
                    if(amb[a][j]==bases[k][0]) {
                        //bComps[k] stores all the rev comps of bases[k]
                        for(int l=0;l<bComps[k].length;l++) {
                            //check the base/code is not used already in cText
                            boolean test=false;
                            for(int q=0;q<cText.length();q++) {
                                if(bComps[k][l]==cText.charAt(q)) {
                                    test=true;
                                }
                            }
                            //if not add it
                            if(!test) {
                                cText+=bComps[k][l];
                                cCount++;
                            }
                        }
                        break;
                    }
                }//end of bases[k] loop
            }//end of amb[i][j] loop

            aComps[a]=new char[cText.length()];
            for(int j=0;j<cText.length();j++) {
                aComps[a][j]=cText.charAt(j);
            }
        }//end of amb[a] loop

        //Basically compsBase stores the actual base or ambiguity code
        //Whilst comps stores all the rev comps [bases and ambiguity] for the compsBase
        comps=new char[bComps.length+aComps.length][];
        for(int i=0;i<bComps.length;i++) {
            comps[i]=bComps[i];
        }
        for(int i=0;i<aComps.length;i++) {
            comps[i+bComps.length]=aComps[i];
        }

        compsBase=new char[comps.length];
        for(int i=0;i<bases.length;i++) {
            compsBase[i]=bases[i][0];
        }
        for(int i=0;i<amb.length;i++) {
            compsBase[i+bases.length]=amb[i][0];
        }


        //Transitions - no need to store transversion - if it is not a Ts it is a Tv
        ts=new char[2][2];
        ts[0][0]='A';ts[0][1]='G';//AG -> purines
        ts[1][0]='C';ts[1][1]='T';//CT -> pyrimidines

        //this is essentially a merger of bases and amb
        //single array to store what each base/amiguity represents
        //the bases represent themselves!
        //used in check mutations for checking TsTv
        muts=new char[bases.length+amb.length][];
        muts[0]=new char[2];muts[1]=new char[2];muts[2]=new char[2];muts[3]=new char[2];
        muts[0][0]='A';muts[0][1]='A';
        muts[1][0]='C';muts[1][1]='C';
        muts[2][0]='G';muts[2][1]='G';
        muts[3][0]='T';muts[3][1]='T';
        for(int i=0;i<amb.length;i++) {
            muts[i+bases.length]=amb[i];//amb[i] is an array
        }

    }//end of setUpBases()


    public  void setUpCts() {

        minPm=1-0.8205;//maximum Primer mismatch % based on combined F and R length
        minPmI=1-0.7;//maximum INDIVIDUAL Primer mismacth % for an individual primer
        minPb=1-0.85;//maximum Probe mismtach %

        max14nt=2;//maximum mutations in 1-4 nts at 3' end of primer


    }//end of setUpCts

    void checkForwardDirection(int s){
        //F/P/R primers in ForwardDirection checked against the ComplementSeqs
        for(int p = 0; p<(primerNum /2); p++) {
            for(int i=0;i<=compSeqs[s].length()-pLens[p];i++) {

                //this is where to record the result - want it on the 3' end of the primer binding site
                int pos=i+pLens[p]-1;

                for(int j=0;j<pLens[p];j++) {

                    if(compSeqs[s].charAt(i+j)=='N')
                        scores[pos][p][14]++;//number of Ns

                    int mut=checkMutation(primers[p].charAt(j),compSeqs[s].charAt(i+j));

                    if(mut>0) {
                        scores[pos][p][0]++;//total mismatches

                        if(mut==1) {
                            scores[pos][p][1]++;//total ts
                        }
                        else if(mut==2) {
                            scores[pos][p][2]++;//total tv
                        }

                        if(j==pLens[p]-1) {
                            if(mut==1)
                                scores[pos][p][3]++;//ts at nt1
                            else if(mut==2)
                                scores[pos][p][4]++;//tv at nt1

                            scores[pos][p][11]++;//mismatches nt1-4
                        }

                        else if(j==pLens[p]-2) {
                            if(mut==1)
                                scores[pos][p][5]++;//ts at nt2
                            else if(mut==2)
                                scores[pos][p][6]++;//tv at nt2

                            scores[pos][p][11]++;//mismatches nt1-4
                        }

                        else if(j==pLens[p]-3) {
                            if(mut==1)
                                scores[pos][p][7]++;//ts at nt3
                            else if(mut==2)
                                scores[pos][p][8]++;//tv at nt3

                            scores[pos][p][11]++;//mismatches nt1-4
                        }

                        else if(j==pLens[p]-4) {
                            if(mut==1)
                                scores[pos][p][9]++;//ts at nt4
                            else if(mut==2)
                                scores[pos][p][10]++;//tv at nt4

                            scores[pos][p][11]++;//mismatches nt1-4
                        }

                        if(j<4) {
                            scores[pos][p][15]++;//mismatches at 5' 1-4nt
                        }

                    }//end of mut>0
                }//end or primer characters loop j
            }//end of sequence characters loop i
        }//end of F/P/R primer (forward) loop p

    }

    void checkReverseDirection(int s){
        //F/P/R primers in ReverseDirection checked against Seqs
        for(int p = (primerNum /2); p<(primerNum); p++) {
            for(int i=0;i<=seqs[s].length()-pLens[p];i++) {

                int pos=i;

                for(int j=0;j<pLens[p];j++) {

                    if(seqs[s].charAt(i+j)=='N')
                        scores[pos][p][14]++;

                    int mut=checkMutation(primers[p].charAt(j),seqs[s].charAt(i+j));

                    if(mut>0) {
                        scores[pos][p][0]++;//total mismatches

                        if(mut==1) {
                            scores[pos][p][1]++;//total ts
                        }
                        else if(mut==2) {
                            scores[pos][p][2]++;//total tv
                        }

                        if(j==0) {
                            if(mut==1)
                                scores[pos][p][3]++;//ts at nt1
                            else if(mut==2)
                                scores[pos][p][4]++;//tv at nt1

                            scores[pos][p][11]++;//mismatches nt1-4
                        }

                        else if(j==1) {
                            if(mut==1)
                                scores[pos][p][5]++;//ts at nt2
                            else if(mut==2)
                                scores[pos][p][6]++;//tv at nt2

                            scores[pos][p][11]++;//mismatches nt1-4
                        }

                        else if(j==2) {
                            if(mut==1)
                                scores[pos][p][7]++;//ts at nt3
                            else if(mut==2)
                                scores[pos][p][8]++;//tv at nt3

                            scores[pos][p][11]++;//mismatches nt1-4
                        }

                        else if(j==3) {
                            if(mut==1)
                                scores[pos][p][9]++;//ts at nt4
                            else if(mut==2)
                                scores[pos][p][10]++;//tv at nt4

                            scores[pos][p][11]++;//mismatches nt1-4
                        }

                        if(j>=pLens[p]-4){
                            scores[pos][p][15]++;//mismactehs at 5' 1-4nt
                        }
                    }//end of mut>0
                }//end or primer characters loop j
            }//end of sequence characters loop i
        }//end of F/P/R primer loop p
    }

    public void setMismatchOverhang(int s){
        //Set the overhangs at the ends to max mismatches - only records score at the 3' position of primer so empty values at the ends
        for(int i=0;i<scores.length;i++) {//i=position
            for(int j=0;j<scores[i].length;j++) {//j=primer
                if(j<(primerNum /2)) {
                    //Matching Forward Direction Primers to 3'-5' stand so empty vals at left 3' end
                    if(i<pLens[j]-1) {
                        scores[i][j][0]=pLens[j];
                    }
                }
                else {
                    //Matching Reverse Direction Primers to 5'-3' stand so empty vals at right 3' end
                    if(i>compSeqs[s].length()-pLens[j]) {
                        scores[i][j][0]=pLens[j];
                    }
                }
            }
        }
    }

    public void ruleOutPos(){

        //Rule out positions based on % mismatch and number of mismatches at 4nt 3' end
        for(int i=0;i<scores.length;i++) {
            for(int j=0;j<scores[i].length;j++) {
                scores[i][j][13] = 0;//0=pass, >0 is fail, set to pass initially

                //check Ns - bit crude - what to discount regions with high N
                double ns = (double) scores[i][j][14] / (double) pLens[j];
                if (ns >= minPb | ns >= minPm)//
                    scores[i][j][13] += 4;

                //probe
                if (hasProbe){
                    if (j == 1 | j == 4) {
                        double perc = (double) scores[i][j][0] / (double) primers[j].length();
                        if (perc >= minPb) {
                            scores[i][j][12] = 1;//flag for failed % test
                            scores[i][j][13] += 2;
                        }
                    }
                }
                //primer
                else {
                    //if more than 2 mismatches in 1-4nt at 3'
                    if(scores[i][j][11]>max14nt) {
                        scores[i][j][13]+=1;
                    }
                    //use mLen (combined F and R length) initially to find initial candidates - filter later
                    double perc=(double)scores[i][j][0]/(double)mLen;
                    double percI=(double)scores[i][j][0]/(double)pLens[j];

                    if(perc>=minPm | percI>=minPmI) {
                        scores[i][j][12]=1;//flag for failed % test
                        scores[i][j][13]+=2;
                    }
                }
            }
        }//end of finding positions that prime loop
    }

    public void PrintCandidatePos(BufferedWriter outOut) throws IOException {
        //count candidate positions
        int can[]=new int[primerNum];

        for(int i=0;i<scores[0].length;i++) {
            count=0;

            for(int j=0;j<scores.length;j++) {
                if(scores[j][i][13]==0) {
                    double mis=100-(double)scores[j][i][0]/(double)pLens[i]*100;
                    double mis2=mis;
                    if(hasProbe) {
                        if (i != 1 & i != 4) {
                            mis2 = 100 - (double) scores[j][i][0] / (double) (mLen) * 100;
                            tTx = (j + 1) + " position is a candidate for " + pNames[i] + " " + df.format(mis) + "% " + df.format(mis2) + "% [%Match %MatchPair], " + scores[j][i][0] + " " + scores[j][i][11] + " [TotMis Tot1-4]";
                            System.out.println(tTx);
                            outOut.write(tTx + System.getProperty("line.separator"));
                        }
                    } else {
                        mis2 = 100 - (double) scores[j][i][0] / (double) (mLen) * 100;
                        tTx = (j + 1) + " position is a candidate for " + pNames[i] + " " + df.format(mis) + "% " + df.format(mis2) + "% [%Match %MatchPair], " + scores[j][i][0] + " " + scores[j][i][11] + " [TotMis Tot1-4]";
                        System.out.println(tTx);
                        outOut.write(tTx + System.getProperty("line.separator"));
                    }
                    count++;
                }
            }
            can[i]=count;
        }
        if(hasProbe) {
            //RJO - at this point the probe is not checked if it is inbetween the two
            //Technically probes could go in any set/direction as long as between the F/R primers
            tTx = can[0] + "-[" + (can[1] + can[4]) + "]-" + can[5] + " Fwd-[Probe]-Rev individual candidate primer/probe positions found in expected orientation";
            System.out.println(tTx);
            outOut.write(tTx + System.getProperty("line.separator"));

            tTx = can[3] + "-[" + (can[1] + can[4]) + "]-" + can[2] + " Fwd-[Probe]-Rev individual candidate primer/probe positions found in opposite orientation";
            System.out.println(tTx);
            outOut.write(tTx + System.getProperty("line.separator"));
        } else{
            //RJO - at this point the probe is not checked if it is inbetween the two
            //Technically probes could go in any set/direction as long as between the F/R primers
            tTx = can[0] + "--" + can[3] + " Fwd-[Probe]-Rev individual candidate primer/probe positions found in expected orientation";
            System.out.println(tTx);
            outOut.write(tTx + System.getProperty("line.separator"));

            tTx = can[2] + "--" + can[1] + " Fwd-[Probe]-Rev individual candidate primer/probe positions found in opposite orientation";
            System.out.println(tTx);
            outOut.write(tTx + System.getProperty("line.separator"));
        }
    }

    public void buildSeqList(EnvMaker maker, int seqCount){

        seqs=new String[seqCount];
        seqNames=new String[seqCount];
        for(int i=0;i<seqs.length;i++) {
            seqs[i]="";
            seqNames[i]="";
        }
        count=0;
        File inFile = new File(sFilename);
        try {
            BufferedReader input =  new BufferedReader(new FileReader(inFile));

            try {
                String line = null;

                while (( line = input.readLine()) != null) {
                    if(line.indexOf(">")==0) {
                        count++;
                        seqNames[count-1]=line;
                    }
                    else {
                        seqs[count-1]+=line;
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

        System.out.println("");
        System.out.println(count+" sequences found in target sequence file");
        for(int i=0;i<seqs.length;i++) {
            seqs[i]= maker.checkSeqs(seqNames[i], seqs[i]);
        }

        compSeqs=new String[seqs.length];
        for(int i=0;i<seqs.length;i++) {
            compSeqs[i]=complement(seqs[i]);
        }

        System.out.println("");
    }

    public void runMismatchDetection(int s){
        checkForwardDirection(s);

        checkReverseDirection(s);

        setMismatchOverhang(s);

        ruleOutPos();

    }
}
