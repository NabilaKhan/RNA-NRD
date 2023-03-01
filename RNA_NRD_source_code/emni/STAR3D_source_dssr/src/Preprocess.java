import java.io.*;
import java.net.URL;
import java.net.HttpURLConnection;
import java.nio.file.*;
import java.util.*;

public class Preprocess {
	public static void main(String[] args) throws Exception {
		if (args.length != 2) {
			System.err.println("usage: java -cp STAR3D.jar Preprocess PDB chain");
			System.exit(0);
		}

		String PDBID = args[0];
		String chainID = args[1];

		int pdb_type = -1;//0 pdb, 1 pdbx

		//get the path STAR3D.jar
		String STAR3D_PATH = new File(Preprocess.class.getProtectionDomain().getCodeSource().getLocation().getPath()).getParentFile().getPath();
		// String STAR3D_PATH = "/home/xiaoli/software/ori_star/ori/STAR3D_source/";//debug

		File PDB_DATA_PATH = new File(STAR3D_PATH, "PDB");
		File SI_DATA_PATH = new File(STAR3D_PATH, "STAR3D_struct_info");

		if (!(PDB_DATA_PATH.exists() && PDB_DATA_PATH.isDirectory())) PDB_DATA_PATH.mkdir();
		if (!(SI_DATA_PATH.exists() && SI_DATA_PATH.isDirectory())) SI_DATA_PATH.mkdir();

		PDBID = PDBID.toLowerCase();
		//chain ID can be lower case
		//chainID=chainID.toUpperCase();

		File PDB_fn = new File(PDB_DATA_PATH, PDBID + ".pdb");
		File PDBx_fn = new File(PDB_DATA_PATH, PDBID + ".cif");
		URL link = null;
		HttpURLConnection connection = null;

		if ((!PDB_fn.exists() || !PDB_fn.isFile()) && (!PDBx_fn.exists() || !PDBx_fn.isFile())) {

			System.out.println("Downloading PDB file " + PDBID + ".pdb...");
			link = new URL("https://files.rcsb.org/download/" + PDBID + ".pdb");
			connection = (HttpURLConnection) link.openConnection();
			connection.setRequestMethod("GET");
			connection.connect();

			int pdb_code = connection.getResponseCode();
			if (pdb_code > 400) {
				System.out.println("Did not find " + PDBID + ".pdb on http://www.rcsb.org. Trying to download .cif file ");

				link = new URL("https://files.rcsb.org/download/" + PDBID + ".cif");
				connection = (HttpURLConnection) link.openConnection();
				connection.setRequestMethod("GET");
				connection.connect();
				int pdbx_code = connection.getResponseCode();
				if (pdbx_code > 400)
					System.out.println("Did not find " + PDBID + ".cif on http://www.rcsb.org.");
				else
					pdb_type = 1;//find pdbx file
			} else
				pdb_type = 0;//find pdb file

			InputStream in = new BufferedInputStream(link.openStream());
			ByteArrayOutputStream out = new ByteArrayOutputStream();
			byte[] buf = new byte[1024];
			int n = 0;
			while (-1 != (n = in.read(buf))) out.write(buf, 0, n);
			out.close();
			in.close();
			byte[] response = out.toByteArray();

			FileOutputStream fos = null;
			if (pdb_type == 0)
				fos = new FileOutputStream(PDB_fn);
			if (pdb_type == 1)
				fos = new FileOutputStream(PDBx_fn);
			fos.write(response);
			fos.close();
		} else {
			if (PDB_fn.exists() && PDB_fn.isFile())
				pdb_type = 0;
			if (PDBx_fn.exists() && PDBx_fn.isFile())
				pdb_type = 1;

		}

		File anno_file = null;
		URL anno_link = null;
		HttpURLConnection anno_connection = null;

		anno_file = new File(SI_DATA_PATH, PDBID + ".dssr");
		if (!anno_file.exists() || anno_file.length() == 0) {
			if (pdb_type == 1) {
				System.out.println("Processing " + PDBID + ".cif with DSSR...");
				PDB_fn = PDBx_fn;
			}
			if (pdb_type == 0)
				System.out.println("Processing " + PDBID + ".pdb with DSSR...");
			

			// ******************* Modified code to download DSSR annotation from website *******************
			try {	
			String pdb_id_low = PDBID.toLowerCase();

			anno_link = new URL("http://skmatic.x3dna.org/pdb/" + pdb_id_low + "/" + pdb_id_low + ".out");
			anno_connection.setFollowRedirects(false);
			anno_connection.setInstanceFollowRedirects(false);
			anno_connection = (HttpURLConnection) anno_link.openConnection();
			anno_connection.setRequestMethod("GET");
			anno_connection.connect();
			int anno_code = anno_connection.getResponseCode();
			System.out.println(anno_code);
			if (anno_code > 400)
				System.out.println("Did not find " + PDBID + ".out on http://skmatic.x3dna.org.");

			InputStream inputStream = new BufferedInputStream(anno_link.openStream());
			Files.copy(inputStream, Paths.get("STAR3D_struct_info/" + pdb_id_low + ".dssr"), StandardCopyOption.REPLACE_EXISTING);


			// InputStream inputStream = new URL("http://skmatic.x3dna.org/pdb/" + pdb_id_low + "/" + pdb_id_low + ".out").openStream();
			// Files.copy(inputStream, Paths.get("STAR3D_struct_info/" + pdb_id_low + ".dssr"), StandardCopyOption.REPLACE_EXISTING);

			} catch (Exception e){
				System.out.println("Missing annotation for PDB file " + PDBID + " on database http://skmatic.x3dna.org.");
                System.out.println("Generate DSSR annotation manually for the PDB file using tool DSSR and put inside folder 'RNA_NRD_source_code/STAR3D_source_dssr/STAR3D_struct_info'. Then try to run the code again!");
                System.exit(0);
			}
				

			//Process p = Runtime.getRuntime().exec(new File(STAR3D_PATH, "tools/DSSR/x3dna-dssr").toString() + " -i=" + PDB_fn);
			//BufferedReader in = new BufferedReader(new InputStreamReader(p.getInputStream()));
			//PrintWriter out = new PrintWriter(anno_file);
			//String line;
			//while ((line = in.readLine()) != null) out.println(line);
			//out.close();
			//in.close();
			//p.waitFor();

			// p = Runtime.getRuntime().exec(new File(STAR3D_PATH, "tools/DSSR/x3dna-dssr").toString() + " --cleanup");
			//p.waitFor();
		}

		File npk_ct_fn = new File(SI_DATA_PATH, PDBID + '_' + chainID + ".npk.ct");
		if (!npk_ct_fn.exists() || npk_ct_fn.length() == 0) {
			System.out.println("Removing pseudoknots from " + PDBID + "_" + chainID + "...");

			PDBParser pdb_parser = null;
			if (pdb_type == 0) {
				pdb_parser = new PDB(PDB_fn, chainID);
			}
			if (pdb_type == 1) {
				pdb_parser = new PDBx(PDBx_fn, chainID);
			}

			String chain_seq = pdb_parser.get_chain_seq(chainID);
			//get the base pairs
			HashSet<Pair<Integer, String>> chain_bp;

			HashSet<Pair<ResID, String>> basepair = new HashSet<Pair<ResID, String>>();

			DSSR dssr_parser = new DSSR(anno_file);

			basepair = dssr_parser.basepair;

			chain_bp = Lib.get_seq_bp(pdb_parser, basepair, chainID);
			Object[][] m = new Object[chain_seq.length()][6];

			for (int i = 0; i < chain_seq.length(); i++) {
				m[i][0] = new Integer(i + 1);
				m[i][1] = new Character(chain_seq.charAt(i));
				m[i][2] = new Integer(i);
				m[i][3] = new Integer(i + 2);
				m[i][4] = new Integer(0);
				m[i][5] = new Integer(i + 1);
			}

			List<Integer> multi_pair = new ArrayList<Integer>();

			//number of Watson-Crick base pairs
			int num_WC = 0;
			for (Pair<Integer, String> P : chain_bp) {
				if (P.v.equals("WWc")) {
					if ((Integer) m[P.left][4] != 0) multi_pair.add(P.left);
					m[P.left][4] = P.right + 1;
					if ((Integer) m[P.right][4] != 0) multi_pair.add(P.right);
					m[P.right][4] = P.left + 1;

					num_WC += 1;
				}
			}

			File ct_fn = new File(SI_DATA_PATH, PDBID + '_' + chainID + ".ct");
			PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(ct_fn)));
			out.printf("%d\t%s_%s\n", chain_seq.length(), PDBID, chainID);
			for (int i = 0; i < chain_seq.length(); i++) {
				out.printf("%d\t%s\t%d\t%d\t%d\t%d\n", m[i][0], m[i][1], m[i][2], m[i][3], m[i][4], m[i][5]);
			}
			out.close();

			//if no pairing interaction, RemovePseudoknot will not work correctly.
			//so we just copy ct to npk.ct, if there is not base pair.
			if (num_WC == 0) {
				Files.copy(ct_fn.toPath(), npk_ct_fn.toPath(), StandardCopyOption.REPLACE_EXISTING);
				return;
			}

			ProcessBuilder pb = new ProcessBuilder();
			Map<String, String> env = pb.environment();
			env.put("DATAPATH", "tools/data_tables");
			Process p = Runtime.getRuntime().exec(new File(STAR3D_PATH, "tools/RemovePseudoknots").toString() + " -m " + ct_fn + " " + npk_ct_fn);
			p.waitFor();

		}
	}
}
