import java.io.*;
import java.util.*;
import java.math.*;


/*
This program reads Gaussian output files and pulls out the NMR shielding constants for each atom in each conformer.  It also extracts energies and calculates a Boltzmann-averaged shielding constant for each atom.

Input syntax is the names of the output files without the .out:

java nmrPredictGaussian myMoleculeCnf001 myMoleculeCnf002 myMoleculeCnf003 ....

(or java nmrPredictGaussian `ls myMoleculeCnf*.out | cut -d"." -f1`)

The program prints the output as a comma separated list suitable for pasting into a spreadsheet.

*/


class nmrPredictGaussian
{

  static double gasConstant = 8.3145;
  static double temperature = 298.15;
  static double hartreeEnergy = 2625.499629554010;


  public static void main(String[] args)
  {
  int nstructs = args.length;
  int natoms = 0;

  String[] fileName = new String[nstructs];
  for(int i=0;i<nstructs;i++)
    {
    fileName[i]=(args[i]+".out");
    }

   BufferedReader buffRead;
   String currLine="";

   try
     {
     buffRead = new BufferedReader(new FileReader(new File(fileName[0])));
     while(buffRead.ready())
       {
       currLine=buffRead.readLine();
       if(currLine.indexOf("Isotropic")!=-1)
         {
         natoms++;
         }
       }
     buffRead.close();
     }
   catch(Exception e)
     {
     e.printStackTrace();
     }
 
    String[] atomLabels = new String[natoms];
    int[] atomNumbers = new int[natoms];
    double[][] shieldingConstants = new double[nstructs][natoms];
    double[] energies = new double[nstructs];
    double[] relenergies = new double[nstructs];
    double[] populations = new double[nstructs];
    double[] boltzshifts = new double[natoms];

    currLine="";

    for(int i=0;i<nstructs;i++)
      {
      try
        {
        buffRead = new BufferedReader(new FileReader(new File(fileName[i])));
        int currAtom=0;
        while(buffRead.ready())
          {
          currLine=buffRead.readLine();
          if(currLine.indexOf("SCF Done")!=-1)
            {
            StringTokenizer st = new StringTokenizer(currLine);
            for(int j=0;j<4;j++)
              {
              String temp = st.nextToken();
              }
            energies[i] = Double.valueOf(st.nextToken());
            }
          if(currLine.indexOf("Isotropic")!=-1)
            {
            StringTokenizer st = new StringTokenizer(currLine);
            if(i==0)
              {
              atomNumbers[currAtom]=Integer.valueOf(st.nextToken());
              atomLabels[currAtom]=st.nextToken();
              }
            else 
              {
              for(int j=0;j<2;j++) 
                {
                String temp = st.nextToken();
                }
              }
            for(int j=0;j<2;j++)
              {
              String temp = st.nextToken();
              }
            shieldingConstants[i][currAtom] = Double.valueOf(st.nextToken());
            currAtom++;
            }
          }
        buffRead.close();
        }
      catch(Exception e)
        {
        e.printStackTrace();
        }

      }


    //Find global minimum energy
    double minEnergy = energies[0];
    for(int i=0;i<nstructs;i++)
      {
      if(energies[i]<minEnergy)
        {
        minEnergy = energies[i];
        }
      }

    //Fill in relenergies
    for(int i=0;i<nstructs;i++)
      {
      relenergies[i]=(energies[i]-minEnergy)*hartreeEnergy;
      }

    //Calculate populations
    for(int i=0;i<nstructs;i++)
      {
      populations[i]=Math.exp(-relenergies[i]*1000/(gasConstant*temperature));
      }
    double q = 0;
    for(int i=0;i<nstructs;i++)
      {
      q=q+populations[i];
      }
    for(int i=0;i<nstructs;i++)
      {
      populations[i]=populations[i]/q*100;
      }

    //Calculate averaged shifts
    for(int i=0;i<natoms;i++)
      {
      boltzshifts[i]=0;
      for(int j=0;j<nstructs;j++)
        {
        boltzshifts[i]=boltzshifts[i]+shieldingConstants[j][i]*populations[j]/100;
        }
      }




     System.out.print("Gaussian results for,");
     for(int i=0;i<nstructs;i++)
       {
       System.out.print(args[i]+",");
       }
     System.out.print("Boltzmann average"+",");
     System.out.print("\n");
     System.out.print("Energy / hartree,");
     for(int i=0;i<nstructs;i++)
       {
       System.out.print(energies[i]+",");
       }
     System.out.print("\n");
     System.out.print("Relative energy / kJ mol-1,");
     for(int i=0;i<nstructs;i++)
        {
        System.out.print(relenergies[i]+",");
        }
     System.out.print("\n");
     System.out.print("Population (%),");
     for(int i=0;i<nstructs;i++)
        {
        System.out.print(populations[i]+",");
        }
     System.out.print("\n");
     for(int j=0;j<natoms;j++)
       {
       //if((atomLabels[j].compareTo("C")==0)||(atomLabels[j].compareTo("H")==0))
         {
         System.out.print("Shielding Constant for "+atomLabels[j]+atomNumbers[j]+",");
         for(int i=0;i<nstructs;i++)
           {
           System.out.print(shieldingConstants[i][j]+",");
           }
         System.out.print(boltzshifts[j]);
         System.out.print("\n");
         }
       }
     }     




  }
