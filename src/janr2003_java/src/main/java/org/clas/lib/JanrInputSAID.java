package org.clas.lib;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

public class JanrInputSAID extends Constants {
	
	String path = janrPath+"multipols/";

	public JanrInputSAID() {
		
	}
	
	public void loadMultipoles() {
        System.out.println("JanrInputSAID: Reading SAID tables");
		loadMultipoleData("mp33",51,6,mp33r,mp33i);
		loadMultipoleData( "s11",51,9,s11r,s11i);
		loadMultipoleData( "s31",51,9,s31r,s31i);
		loadMultipoleData( "p11",51,9,p11r,p11i);
		loadMultipoleData( "p13",51,9,p13r,p13i);
		loadMultipoleData( "p31",51,9,p31r,p31i);
		loadMultipoleData( "p33",51,9,p33r,p33i);
		loadMultipoleData( "d13",51,9,d13r,d13i);
		loadMultipoleData( "d15",51,9,d15r,d15i);
		loadMultipoleData( "d33",51,9,d33r,d33i);
		loadMultipoleData( "f15",51,9,f15r,f15i);
		loadMultipoleData( "f37",51,9,f37r,f37i);
		loadMultipoleData( "pim",13,3,pimr,pimi);
		loadMultipoleData( "sim",13,3,simr,simi);
	}

	public void loadMultipoleData(String mpol, int nmax, int len, double re[], double im[]) {
        try (BufferedReader br = new BufferedReader(new FileReader(path+mpol+".dat"))) {
        	String line; int i=0; 
            while(i < nmax) {
                line = br.readLine();
                if (line.trim().length() ==  0) continue;
                String[] tokens = line.trim().split("\\s+");
                if (tokens.length < len) break;
                wsa[i]= Double.parseDouble(tokens[0]) / 1000.0;
                re[i] = Double.parseDouble(tokens[1]) / (mpol=="mp33"?52.437:1.0);
                im[i] = Double.parseDouble(tokens[2]) / (mpol=="mp33"?52.437:1.0); i++;
                // at2, at3, at4 are not used
            }
        } catch (IOException e) {
            e.printStackTrace();
        }		
	}

    public static void main(String[] args) {
        JanrInputSAID loader = new JanrInputSAID();
        loader.loadMultipoles();
    }
}
