package org.clas.lib;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

public class JanrInputSAID extends Constants {

	public JanrInputSAID() {
		
	}

    public void loadMultipoleData()  {
        System.out.println("JANR_INPUT_SAID: Reading SAID tables");

        // mp33.dat
        try (BufferedReader br = new BufferedReader(new FileReader("multipols/mp33.dat"))) {
            for (int i = 0; i < N51; i++) {
                String line = br.readLine();
                if (line == null) break;
                String[] tokens = line.trim().split("\\s+");
                if (tokens.length < 6) break;
                wsa[i]   = Double.parseDouble(tokens[0]) / 1000.0;
                mp33r[i] = Double.parseDouble(tokens[1]) / 52.437;
                mp33i[i] = Double.parseDouble(tokens[2]) / 52.437;
                // at2, at3, at4 are not used
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        // s11.dat
        try (BufferedReader br = new BufferedReader(new FileReader("multipols/s11.dat"))) {
            for (int i = 0; i < N51; i++) {
                String line = br.readLine();
                if (line == null) break;
                String[] tokens = line.trim().split("\\s+");
                if (tokens.length < 9) break;
                wsa[i]  = Double.parseDouble(tokens[0]) / 1000.0;
                s11r[i] = Double.parseDouble(tokens[5]);
                s11i[i] = Double.parseDouble(tokens[6]);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        // s31.dat
        try (BufferedReader br = new BufferedReader(new FileReader("multipols/s31.dat"))) {
            for (int i = 0; i < N51; i++) {
                String line = br.readLine();
                if (line == null) break;
                String[] tokens = line.trim().split("\\s+");
                if (tokens.length < 9) break;
                wsa[i]  = Double.parseDouble(tokens[0]) / 1000.0;
                s31r[i] = Double.parseDouble(tokens[5]);
                s31i[i] = Double.parseDouble(tokens[6]);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        // p11.dat
        try (BufferedReader br = new BufferedReader(new FileReader("multipols/p11.dat"))) {
            for (int i = 0; i < N51; i++) {
                String line = br.readLine();
                if (line == null) break;
                String[] tokens = line.trim().split("\\s+");
                if (tokens.length < 9) break;
                wsa[i]  = Double.parseDouble(tokens[0]) / 1000.0;
                p11r[i] = Double.parseDouble(tokens[5]);
                p11i[i] = Double.parseDouble(tokens[6]);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        // p13.dat
        try (BufferedReader br = new BufferedReader(new FileReader("multipols/p13.dat"))) {
            for (int i = 0; i < N51; i++) {
                String line = br.readLine();
                if (line == null) break;
                String[] tokens = line.trim().split("\\s+");
                if (tokens.length < 9) break;
                wsa[i]  = Double.parseDouble(tokens[0]) / 1000.0;
                p13r[i] = Double.parseDouble(tokens[5]);
                p13i[i] = Double.parseDouble(tokens[6]);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        // p31.dat
        try (BufferedReader br = new BufferedReader(new FileReader("multipols/p31.dat"))) {
            for (int i = 0; i < N51; i++) {
                String line = br.readLine();
                if (line == null) break;
                String[] tokens = line.trim().split("\\s+");
                if (tokens.length < 9) break;
                wsa[i]  = Double.parseDouble(tokens[0]) / 1000.0;
                p31r[i] = Double.parseDouble(tokens[5]);
                p31i[i] = Double.parseDouble(tokens[6]);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        // p33.dat
        try (BufferedReader br = new BufferedReader(new FileReader("multipols/p33.dat"))) {
            for (int i = 0; i < N51; i++) {
                String line = br.readLine();
                if (line == null) break;
                String[] tokens = line.trim().split("\\s+");
                if (tokens.length < 9) break;
                wsa[i]  = Double.parseDouble(tokens[0]) / 1000.0;
                p33r[i] = Double.parseDouble(tokens[5]);
                p33i[i] = Double.parseDouble(tokens[6]);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        // d13.dat
        try (BufferedReader br = new BufferedReader(new FileReader("multipols/d13.dat"))) {
            for (int i = 0; i < N51; i++) {
                String line = br.readLine();
                if (line == null) break;
                String[] tokens = line.trim().split("\\s+");
                if (tokens.length < 9) break;
                wsa[i]  = Double.parseDouble(tokens[0]) / 1000.0;
                d13r[i] = Double.parseDouble(tokens[5]);
                d13i[i] = Double.parseDouble(tokens[6]);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        // d15.dat
        try (BufferedReader br = new BufferedReader(new FileReader("multipols/d15.dat"))) {
            for (int i = 0; i < N51; i++) {
                String line = br.readLine();
                if (line == null) break;
                String[] tokens = line.trim().split("\\s+");
                if (tokens.length < 9) break;
                wsa[i]  = Double.parseDouble(tokens[0]) / 1000.0;
                d15r[i] = Double.parseDouble(tokens[5]);
                d15i[i] = Double.parseDouble(tokens[6]);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        // d33.dat
        try (BufferedReader br = new BufferedReader(new FileReader("multipols/d33.dat"))) {
            for (int i = 0; i < N51; i++) {
                String line = br.readLine();
                if (line == null) break;
                String[] tokens = line.trim().split("\\s+");
                if (tokens.length < 9) break;
                wsa[i]  = Double.parseDouble(tokens[0]) / 1000.0;
                d33r[i] = Double.parseDouble(tokens[5]);
                d33i[i] = Double.parseDouble(tokens[6]);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        // f15.dat
        try (BufferedReader br = new BufferedReader(new FileReader("multipols/f15.dat"))) {
            for (int i = 0; i < N51; i++) {
                String line = br.readLine();
                if (line == null) break;
                String[] tokens = line.trim().split("\\s+");
                if (tokens.length < 9) break;
                wsa[i]  = Double.parseDouble(tokens[0]) / 1000.0;
                f15r[i] = Double.parseDouble(tokens[5]);
                f15i[i] = Double.parseDouble(tokens[6]);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        // f37.dat
        try (BufferedReader br = new BufferedReader(new FileReader("multipols/f37.dat"))) {
            for (int i = 0; i < N51; i++) {
                String line = br.readLine();
                if (line == null) break;
                String[] tokens = line.trim().split("\\s+");
                if (tokens.length < 9) break;
                wsa[i]  = Double.parseDouble(tokens[0]) / 1000.0;
                f37r[i] = Double.parseDouble(tokens[5]);
                f37i[i] = Double.parseDouble(tokens[6]);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        // pim.dat
        try (BufferedReader br = new BufferedReader(new FileReader("multipols/pim.dat"))) {
            for (int i = 0; i < N13; i++) {
                String line = br.readLine();
                if (line == null) break;
                String[] tokens = line.trim().split("\\s+");
                if (tokens.length < 3) break;
                wsa[i]  = Double.parseDouble(tokens[0]) / 1000.0;
                pimr[i] = Double.parseDouble(tokens[1]);
                pimi[i] = Double.parseDouble(tokens[2]);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        // sim.dat
        try (BufferedReader br = new BufferedReader(new FileReader("multipols/sim.dat"))) {
            for (int i = 0; i < N13; i++) {
                String line = br.readLine();
                if (line == null) break;
                String[] tokens = line.trim().split("\\s+");
                if (tokens.length < 3) break;
                wsa[i]  = Double.parseDouble(tokens[0]) / 1000.0;
                simr[i] = Double.parseDouble(tokens[1]);
                simi[i] = Double.parseDouble(tokens[2]);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    // For testing
    public static void main(String[] args) {
        JanrInputSAID loader = new JanrInputSAID();
        loader.loadMultipoleData();
    }
}
