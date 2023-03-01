import java.io.*;
import java.util.*;
import java.util.concurrent.*;

import java.text.DateFormat;
import java.text.SimpleDateFormat;

import org.apache.commons.cli.*;
import org.ejml.data.*;
import org.ejml.ops.CommonOps;

public class STAR3D{
/*
 * The global parameters.	
 */
	public static String output_fn="output.aln";
	public static double rmsd_cutoff=4.0;
	public static int min_stack_size=3;
	public static double gap_open_penalty=-5.;
	public static double gap_extend_penalty=-2.;
	public static int thread_num=1;
	public static double match_score=3.;
	public static double mismatch_score=0.;
	public static boolean pdb_output=false;

	public static ArrayList<ResID> ResID1_list = new ArrayList<ResID>();
	public static ArrayList<ResID> ResID2_list = new ArrayList<ResID>();
	
	public static void main(String[] args)  throws Exception{
		Options options=new Options();
		
		options.addOption("o", true, "Output alignment file");
		options.addOption("r", true, "RMSD cutoff");
		options.addOption("s", true, "Minimum stack size");
		options.addOption("g", true, "Gap open penalty");
		options.addOption("e", true, "Gap extension penalty");
		options.addOption("t", true, "Number of threads");
		options.addOption("h", false, "help information");
		options.addOption("m", true, "Match score");
		options.addOption("i", true, "Mismatch score");
		options.addOption("p", false, "Output alignment PDB or not");
		
		CommandLineParser parser=new BasicParser();
		CommandLine cmd=null;
		try{
			cmd=parser.parse(options, args);
		}catch(ParseException e){
			System.err.println(e.getMessage());
			System.exit(0);
		}
		
		HelpFormatter formatter = new HelpFormatter();
		
		if(cmd.hasOption("h")==true){
			formatter.printHelp( "java -jar STAR3D.jar <options> PDB1 chain1 PDB2 chain2", options );
			System.exit(0);
		}
		

		if(cmd.getArgList().size()!=4){
			formatter.printHelp( "java -jar STAR3D.jar <options> PDB1 chain1 PDB2 chain2", options );
			System.exit(0);
		}
		
		String PDBID1=cmd.getArgList().get(0).toString().toLowerCase();
		//String chainID1=cmd.getArgList().get(1).toString().toUpperCase();
		String chainID1=cmd.getArgList().get(1).toString();
		String PDBID2=cmd.getArgList().get(2).toString().toLowerCase();
		//String chainID2=cmd.getArgList().get(3).toString().toUpperCase();
		String chainID2=cmd.getArgList().get(3).toString();


		if(cmd.getOptionValue("o")!=null) STAR3D.output_fn=cmd.getOptionValue("o");
		if(cmd.getOptionValue("r")!=null) STAR3D.rmsd_cutoff=Double.parseDouble(cmd.getOptionValue("r"));
		if(cmd.getOptionValue("s")!=null) STAR3D.min_stack_size=Integer.parseInt(cmd.getOptionValue("s"));
		if(cmd.getOptionValue("g")!=null) STAR3D.gap_open_penalty=Double.parseDouble(cmd.getOptionValue("g"));
		if(cmd.getOptionValue("e")!=null) STAR3D.gap_extend_penalty=Double.parseDouble(cmd.getOptionValue("e"));
		if(cmd.getOptionValue("t")!=null) STAR3D.thread_num=Integer.parseInt(cmd.getOptionValue("t"));
		if(cmd.getOptionValue("m")!=null) STAR3D.match_score=Double.parseDouble(cmd.getOptionValue("m"));
		if(cmd.getOptionValue("i")!=null) STAR3D.mismatch_score=Double.parseDouble(cmd.getOptionValue("i"));
		if(cmd.hasOption("p")==true) STAR3D.pdb_output=true;
			
/*
 * Parse MCA files		
 */
		String STAR3D_PATH = new File(STAR3D.class.getProtectionDomain().getCodeSource().getLocation().getPath()).getParentFile().getPath();
		// String STAR3D_PATH = "/home/xiaoli/software/ori_star/ori/STAR3D_source/";//debug

		File PDB_DATA_PATH=new File(STAR3D_PATH, "PDB"); 
		File SI_DATA_PATH=new File(STAR3D_PATH, "STAR3D_struct_info");

		File PDB1_fn = new File(PDB_DATA_PATH, PDBID1 + ".pdb");
		File PDB2_fn = new File(PDB_DATA_PATH, PDBID2 + ".pdb");

		File PDBx1_fn = new File(PDB_DATA_PATH, PDBID1 + ".cif");
		File PDBx2_fn = new File(PDB_DATA_PATH, PDBID2 + ".cif");

		int pdb1_type = -1, pdb2_type = -1;

		if (PDB1_fn.exists() && PDB1_fn.isFile()) pdb1_type = 0;
		else if (PDBx1_fn.exists() && PDBx1_fn.isFile()) pdb1_type = 1;

		if (PDB2_fn.exists() && PDB2_fn.isFile()) pdb2_type = 0;
		else if (PDBx2_fn.exists() && PDBx2_fn.isFile()) pdb2_type = 1;

		PDBParser PDB1_parser = null;
		if (pdb1_type == 0)
			PDB1_parser = new PDB(PDB1_fn, chainID1);
		else
			PDB1_parser = new PDBx(PDBx1_fn, chainID1);

		PDBParser PDB2_parser = null;
		if (pdb2_type == 0)
			PDB2_parser = new PDB(PDB2_fn, chainID2);
		else
			PDB2_parser = new PDBx(PDBx2_fn, chainID2);

		ArrayList<Residue> Res1_list = PDB1_parser.get_chain_res().get(chainID1);
		ArrayList<Residue> Res2_list = PDB2_parser.get_chain_res().get(chainID2);

		for(Residue R: Res1_list) ResID1_list.add(R.rid);
		for(Residue R: Res2_list) ResID2_list.add(R.rid);

		String seq1 = PDB1_parser.get_chain_seq(chainID1);
		String seq2 = PDB2_parser.get_chain_seq(chainID2);

		ArrayList<Point> coord1 = PDB1_parser.get_chain_centroid(chainID1);
		ArrayList<Point> coord2 = PDB2_parser.get_chain_centroid(chainID2);

//parse the dssr annotation files
		File anno1_file = new File(SI_DATA_PATH, PDBID1 + ".dssr");
		File anno2_file = new File(SI_DATA_PATH, PDBID2 + ".dssr");
//get all information from the annotation files
		DSSR DSSR1_parser = new DSSR(anno1_file);
		DSSR DSSR2_parser = new DSSR(anno2_file);
//get pair info for a specific chain (chainID)		
		HashSet<Pair<Integer, String>> bp1=Lib.get_seq_bp(PDB1_parser, DSSR1_parser.basepair, chainID1);
		HashSet<Pair<Integer, String>> bp2=Lib.get_seq_bp(PDB2_parser, DSSR2_parser.basepair, chainID2);
		HashMap<Integer, ArrayList<Integer>> bp1_map=Lib.get_seq_pairing_map(bp1);
		HashMap<Integer, ArrayList<Integer>> bp2_map=Lib.get_seq_pairing_map(bp2);
/*
 * get base pairs in the npk secondary structure
 */
		File ct1_fn=new File(SI_DATA_PATH, PDBID1+"_"+chainID1+".npk.ct");
		File ct2_fn=new File(SI_DATA_PATH, PDBID2+"_"+chainID2+".npk.ct");
		
		List<Pair<Integer, Integer>> npk_bp1=Lib.get_ct_bp(ct1_fn);
		List<Pair<Integer, Integer>> npk_bp2=Lib.get_ct_bp(ct2_fn);
		
/*
 * use multiple threading to find ungapped stack mapping		
 */
		ExecutorService exec=Executors.newFixedThreadPool(thread_num);
		List<Future<List<Stackmap>>> Stackmap_multi_thread_ret=new ArrayList<Future<List<Stackmap>>>();
	
		for(int i=0; i<thread_num; i++)
			Stackmap_multi_thread_ret.add(exec.submit(new FindStackmap(i, thread_num, npk_bp1, npk_bp2, coord1, coord2)));
		
		List<Stackmap> SM_all=new ArrayList<Stackmap>();
		for(Future<List<Stackmap>> fs: Stackmap_multi_thread_ret){
			try{
				SM_all.addAll(fs.get());
			}catch(InterruptedException e){
				System.out.println(e);
			}catch(ExecutionException e){
				System.out.println(e);
			}finally{
				exec.shutdown();
			}
		}

		if(SM_all.size()==0){
			System.out.println("No similar stacks");
			System.exit(0);
		}
		
		Collections.sort(SM_all);
		List<Stackmap> SM_top=SM_all.subList(0, Math.min(SM_all.size(), 200));
		
/*
 * find stack map configuration
 */		
		exec=Executors.newFixedThreadPool(thread_num);
		List<Future<List<Pair<Integer, Integer>>>> SMC_multi_thread_ret=new ArrayList<Future<List<Pair<Integer, Integer>>>>();
		
		for(int i=0; i<thread_num; i++)
			SMC_multi_thread_ret.add(exec.submit(new FindCompatibleSM(i, thread_num, SM_top, coord1, coord2)));
		
		List<Pair<Integer, Integer>> SMC_all=new ArrayList<Pair<Integer, Integer>>();
		for(Future<List<Pair<Integer, Integer>>> fs: SMC_multi_thread_ret){
			try{
				SMC_all.addAll(fs.get());
			}catch(InterruptedException e){
				System.out.println(e);
			}catch(ExecutionException e){
				System.out.println(e);
			}finally{
				exec.shutdown();
			}
		}
		
		int[][] SMC_graph=new int[SM_top.size()][SM_top.size()];
		for(Pair<Integer, Integer> P: SMC_all) {
			if(P.v==1) {
				SMC_graph[P.left][P.right]=1;
				SMC_graph[P.right][P.left]=1;
			}
		}
		
		Set<Integer> SMC_graph_vertex=new HashSet<Integer>();
		for(int i=0; i<SMC_graph.length; i++)
			SMC_graph_vertex.add(i);
		
		List<Set<Integer>> SM_clique=new ArrayList<Set<Integer>>();
		Lib.Bron_Kerbosch(SMC_graph, new HashSet<Integer>(), SMC_graph_vertex, new HashSet<Integer>(), SM_clique);
		
		int max_clique_size=0;
		List<Set<Integer>> max_cliques=new ArrayList<Set<Integer>>();
		for(Set<Integer> C: SM_clique){
			int BP_size=0;
			for(Integer I: C){
				BP_size+=SM_top.get(I).size;
			}
			if(BP_size==max_clique_size) max_cliques.add(C);
			if(BP_size>max_clique_size){
				max_clique_size=BP_size;
				max_cliques.clear();
				max_cliques.add(C);
			}
		}

/*
 * find the mapping loops between two structures
 */
		
		List<Integer> optimal_Map1_index=new ArrayList<Integer>(); 
		List<Integer> optimal_Map2_index=new ArrayList<Integer>();
		Point optimal_Map1_XC=new Point(0, 0, 0);
		Point optimal_Map2_YC=new Point(0, 0, 0);
		DenseMatrix64F optimal_Map_R=new DenseMatrix64F();
		double optimal_score=0.;
		double optimal_rmsd=256.0;
		//try all the clique with max size
		for(Set<Integer> C: max_cliques){
			List<Pair<Integer, Integer>> Map1_stack=new ArrayList<Pair<Integer, Integer>>();
			List<Pair<Integer, Integer>> Map2_stack=new ArrayList<Pair<Integer, Integer>>();
			//get the stack mapping: Map1_stack, Map2_stack
			for(Integer I: C){
				 	Stackmap SM=SM_top.get(I);
					Map1_stack.add(new Pair<Integer, Integer>(SM.i1, SM.j1, SM.size));
					Map2_stack.add(new Pair<Integer, Integer>(SM.i2, SM.j2, SM.size));
			}
			Collections.sort(Map1_stack);
			Collections.sort(Map2_stack);
			//get the stack tree	
			Snode root1=Snode.stack_tree(Map1_stack, seq1.length());
			Snode root2=Snode.stack_tree(Map2_stack, seq2.length());
			
			List<Pair<Integer, Integer>> Map1_loop=new ArrayList<Pair<Integer, Integer>>();
			List<Pair<Integer, Integer>> Map2_loop=new ArrayList<Pair<Integer, Integer>>();
			//get the loop mapping from the tree
			Snode.order_loop(root1, Map1_loop);
			Snode.order_loop(root2, Map2_loop);
			
			//get the mapping between stack residues
			List<Integer> Map1_stack_index=new ArrayList<Integer>();
			List<Integer> Map2_stack_index=new ArrayList<Integer>();
			for(Pair<Integer, Integer> P: Map1_stack){
				for(int i=P.left; i<P.left+P.v; i++) Map1_stack_index.add(i);
				for(int i=P.right; i>P.right-P.v; i--) Map1_stack_index.add(i);
			}
			for(Pair<Integer, Integer> P: Map2_stack){
				for(int i=P.left; i<P.left+P.v; i++) Map2_stack_index.add(i);
				for(int i=P.right; i>P.right-P.v; i--) Map2_stack_index.add(i);
			}

			Collections.sort(Map1_stack_index);
			Collections.sort(Map2_stack_index);
			
			//get the transition and rotate matrix for the stacking mapping
			DenseMatrix64F Map1_stack_X=new DenseMatrix64F(Map1_stack_index.size(), 3);
			DenseMatrix64F Map2_stack_Y=new DenseMatrix64F(Map2_stack_index.size(), 3);
			
			for(int i=0; i<Map1_stack_index.size(); i++){
				Map1_stack_X.set(i, 0, coord1.get(Map1_stack_index.get(i)).x);
				Map1_stack_X.set(i, 1, coord1.get(Map1_stack_index.get(i)).y);
				Map1_stack_X.set(i, 2, coord1.get(Map1_stack_index.get(i)).z);
				
				Map2_stack_Y.set(i, 0, coord2.get(Map2_stack_index.get(i)).x);
				Map2_stack_Y.set(i, 1, coord2.get(Map2_stack_index.get(i)).y);
				Map2_stack_Y.set(i, 2, coord2.get(Map2_stack_index.get(i)).z);
			}
			
			Point stack_XC=Geom.centroid(Map1_stack_X);
			Point stack_YC=Geom.centroid(Map2_stack_Y);
			
			DenseMatrix64F stack_R=new DenseMatrix64F();
			stack_R=Geom.Kabsch(Geom.translation(Map1_stack_X, stack_XC), Geom.translation(Map2_stack_Y, stack_YC));	//stack rotation matrix
			
			//get the mapping between the loop residues
			exec=Executors.newFixedThreadPool(thread_num);
			List<Future<List<Pair<Integer, Integer>>>> Align_multi_thread_ret=new ArrayList<Future<List<Pair<Integer, Integer>>>>();
			
			for(int i=0; i<thread_num; i++)
				Align_multi_thread_ret.add(exec.submit(new LoopAlign(i, thread_num, Map1_loop, Map2_loop, bp1_map, bp2_map, coord1, coord2, stack_R, stack_XC, stack_YC)));
			
			List<Pair<Integer, Integer>> Map_loop_align=new ArrayList<Pair<Integer, Integer>>();
			for(Future<List<Pair<Integer, Integer>>> fs: Align_multi_thread_ret){
				try{
					Map_loop_align.addAll(fs.get());
				}catch(InterruptedException e){
					System.out.println(e);
				}catch(ExecutionException e){
					System.out.println(e);
				}finally{
					exec.shutdown();
				}
			}
			List<Integer> Map1_loop_index=new ArrayList<Integer>();
			List<Integer> Map2_loop_index=new ArrayList<Integer>();
			for(Pair<Integer, Integer> P: Map_loop_align){
				if(P.left!=-1 && P.right!=-1){
					Map1_loop_index.add(P.left);
					Map2_loop_index.add(P.right);
				}
			}
			
			//combine the mapped residues
			List<Integer> Map1_index=new ArrayList<Integer>(Map1_stack_index);
			Map1_index.addAll(Map1_loop_index);
			List<Integer> Map2_index=new ArrayList<Integer>(Map2_stack_index);
			Map2_index.addAll(Map2_loop_index);
			
			Collections.sort(Map1_index);
			Collections.sort(Map2_index);
			
			//get the transition and rotate matrix for all mapped residues
			DenseMatrix64F Map1_X=new DenseMatrix64F(Map1_index.size(), 3);
			DenseMatrix64F Map2_Y=new DenseMatrix64F(Map2_index.size(), 3);
			
			for(int i=0; i<Map1_index.size(); i++){
				Map1_X.set(i, 0, coord1.get(Map1_index.get(i)).x);
				Map1_X.set(i, 1, coord1.get(Map1_index.get(i)).y);
				Map1_X.set(i, 2, coord1.get(Map1_index.get(i)).z);
				
				Map2_Y.set(i, 0, coord2.get(Map2_index.get(i)).x);
				Map2_Y.set(i, 1, coord2.get(Map2_index.get(i)).y);
				Map2_Y.set(i, 2, coord2.get(Map2_index.get(i)).z);
			}
			double Map_rmsd=Geom.superimpose(Map1_X, Map2_Y);
			
			Point Map1_XC=Geom.centroid(Map1_X);	//centroid for mapped residues in X
			Point Map2_YC=Geom.centroid(Map2_Y);	//centroid for mapped residues in Y
			DenseMatrix64F Map1_XM=Geom.translation(Map1_X, Map1_XC);
			DenseMatrix64F Map2_YM=Geom.translation(Map2_Y, Map2_YC);
			
			DenseMatrix64F Map1_XT=new DenseMatrix64F(Map1_XM.numCols, Map1_XM.numRows);
			CommonOps.transpose(Map1_XM, Map1_XT);
			
			DenseMatrix64F Map_R=Geom.Kabsch(Map1_XM, Map2_YM);	//rotate matrix for all mapped residues
			
			DenseMatrix64F Map1_RX=new DenseMatrix64F(3, Map1_XT.numCols);
			CommonOps.mult(Map_R, Map1_XT, Map1_RX);
			CommonOps.transpose(Map1_RX);
			
			int Map_psi=0;
			for(int i=0; i<Map1_index.size(); i++){
				if(Geom.distance(Map1_RX.get(i, 0), Map1_RX.get(i, 1), Map1_RX.get(i, 2), Map2_YM.get(i, 0), Map2_YM.get(i, 1), Map2_YM.get(i, 2))<STAR3D.rmsd_cutoff)
					Map_psi+=1;
			}
			if(optimal_score>-1.*Map_psi+Map_rmsd){
				optimal_score=-1*Map_psi+Map_rmsd;
				optimal_rmsd=Map_rmsd;
				optimal_Map1_index=Map1_index;
				optimal_Map2_index=Map2_index;
				
				optimal_Map1_XC=Map1_XC;
				optimal_Map2_YC=Map2_YC;
				optimal_Map_R=Map_R;
			}
		}
		PrintWriter out=new PrintWriter(new BufferedWriter(new FileWriter(STAR3D.output_fn)));
		DateFormat dateFormat=new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
		Calendar cal=Calendar.getInstance();
		out.println(String.format("#STAR3d alignment for %s and %s (%s)", PDBID1+"_"+chainID1, PDBID2+"_"+chainID2, dateFormat.format(cal.getTime())));
		out.println("############");
		out.println("#Parameters#");
		out.println("############");
		out.println(String.format("#RMSD cutoff: %.1fA", rmsd_cutoff));
		out.println(String.format("#Minimum stack size: %d", min_stack_size));
		out.println(String.format("#Gap open penalty: %.1f", gap_open_penalty));
		out.println(String.format("#Gap extension penalty: %.1f", gap_extend_penalty));
		out.println(String.format("#Match score: %.1f", match_score));
		out.println(String.format("#Mismatch score: %.1f", mismatch_score));
		out.println("#########");
		out.println("#Results#");
		out.println("#########");
		out.println(String.format("#Aligned nucleotide: %d", optimal_Map1_index.size()));
		out.println(String.format("#Alignment RMSD: %.2fA", optimal_rmsd));
		out.println("#Nucleotide mapping:");
		for(int i=0; i<optimal_Map1_index.size(); i++){
			out.println(ResID1_list.get(optimal_Map1_index.get(i))+"<->"+ResID2_list.get(optimal_Map2_index.get(i)));
		}
		out.close();
		
		if(STAR3D.pdb_output==true){
			DenseMatrix64F Atom_coord=new DenseMatrix64F(1, 3);
			DenseMatrix64F Atom_coord_T=new DenseMatrix64F(3, 1);
			DenseMatrix64F Atom_coord_R=new DenseMatrix64F(3, 1);
			out=new PrintWriter(new BufferedWriter(new FileWriter(STAR3D.output_fn+".pdb")));
		
			out.println("MODEL        1");
			for(Atom A: PDB1_parser.chain_atom.get(chainID1.charAt(0))){
				Atom_coord.set(0, 0, A.coord.x);
				Atom_coord.set(0, 1, A.coord.y);
				Atom_coord.set(0, 2, A.coord.z);
				Atom_coord=Geom.translation(Atom_coord, optimal_Map1_XC);
				
				CommonOps.transpose(Atom_coord, Atom_coord_T);
				CommonOps.mult(optimal_Map_R, Atom_coord_T, Atom_coord_R);
				
				out.println("ATOM  "+String.format("%5d", A.sn)+" "+String.format("%-4s", A.atom)+" "+String.format("%3s", A.res.symbol)+" "+
							"A"+String.format("%4d", A.res.rid.seqnum)+A.res.rid.icode+"   "+
							String.format("%8.3f", Atom_coord_R.get(0,0))+String.format("%8.3f", Atom_coord_R.get(1,0))+String.format("%8.3f", Atom_coord_R.get(2,0))+"  1.00 99.99"
							);
			}
			out.println("ENDMDL");
		
			out.println("MODEL        2");
			for(Atom A: PDB2_parser.chain_atom.get(chainID2.charAt(0))){
				Atom_coord.set(0, 0, A.coord.x);
				Atom_coord.set(0, 1, A.coord.y);
				Atom_coord.set(0, 2, A.coord.z);
				Atom_coord=Geom.translation(Atom_coord, optimal_Map2_YC);
				
				out.println("ATOM  "+String.format("%5d", A.sn)+" "+String.format("%-4s", A.atom)+" "+String.format("%3s", A.res.symbol)+" "+
							"B"+String.format("%4d", A.res.rid.seqnum)+A.res.rid.icode+"   "+
							String.format("%8.3f", Atom_coord.get(0,0))+String.format("%8.3f", Atom_coord.get(0,1))+String.format("%8.3f", Atom_coord.get(0,2))+"  1.00 99.99"
							);
			}
			
			out.println("ENDMDL");
		
			out.close();
		}
	}
}
