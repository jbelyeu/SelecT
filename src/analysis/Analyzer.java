package analysis;


import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.NavigableSet;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import errors.FileParsingException;
import errors.IllegalInputException;
import calc.HaplotypeTests;
import tools.*;

public class Analyzer {

	/**
	 * Gets all the data from all the stats output files
	 * If directly called from main it will print the data (any combination of tests) into one big file
	 * If called from Analyzer it will return List<WindowStats> for further analysis
	 * 
	 * FOR WINDOW COMBINING AND NO STATS ANALYSIS
	 * @param args0		window combining flag (
	 * @param args0		out dir
	 * @param args1		stats string ("i:h:x:d:f" or any combination of that)
	 * @param args2		window range ("x-y")
	 * 
	 * FOR STATS ANALYSIS ONLY
	 * @param args0		out dir
	 * @param args1		simulation dir (must have files neutral_simulation.tsv and selection_simulation.tsv)
	 */
	public static void main(String[] args) {
		
		Log log = new Log(Log.type.analysis);
		System.out.println("Starting the final phase of SelecT");
		
		String type = args[0];
		
		try {
			
			if(type.equals("combine")) {
				
				checkCombineArgs(args, log);
				 
				Combiner c = new Combiner(args[1], log);
				c.combineWindows(args[2], args[3]);
				c.writeStats();
			}
			else {
				
				checkStatsAnalArgs(args, log);
				
				Combiner c = new Combiner(args[0], log);
				c.combineWindows();
				
				Analyzer a = new Analyzer(log, new File(args[0]), new File(args[1]), c.getAllStats());
				a.runCmsAnalysis();
				a.writeStats();
			}
			 
		} catch (Exception e) {
			
			System.out.println("CMS Died Prematurely." 
					+ " Check log output for troubleshooting.");
			
			log.addLine("\n\nCMS Died Prematurely. Error in computation.");
			
			e.printStackTrace();
		}

	}
	
	private static void checkStatsAnalArgs(String[] args, Log log) 
			throws IllegalInputException {
		
		if(args.length != 2) {
			String msg = "Error: Parameter length incorrect";
			throw new IllegalInputException(log, msg);
		}
		
		log.addLine("Working Parameters");
		log.addLine("Output Dir:\t" + args[0]);
		log.addLine("Sims Dir:\t" + args[1]);
	}
	
	private static void checkCombineArgs(String[] args, Log log) 
			throws IllegalInputException {
		
		if(args.length != 4) {
			String msg = "Error: Parameter length incorrect";
			throw new IllegalInputException(log, msg);
		}
		
		log.addLine("Working Parameters");
		log.addLine("Output Dir:\t" + args[1]);
		log.addLine("Stats String:\t" + args[2]);
		log.addLine("Window Range:\t" + args[3] + "\n");
	}
	
	
	private final int NUM_TESTS = 5;
	private final double DAF_CUTOFF = 0.2;
	
	private File sim_dir;
	private File out_dir;
	
	private List<WindowStats> all_ws;
	private TreeMap<SNP, Double> cms_scores_prod;//product of scores
	private TreeMap<SNP, Double> cms_scores_mean;//mean of scores
	
	Log log;
	
	public Analyzer(Log log, File out_dir, File sim_dir, List<WindowStats> all_ws) {
		
		this.all_ws = all_ws;
		this.out_dir = out_dir;
		this.sim_dir = sim_dir;
		
		cms_scores_prod = new TreeMap<SNP, Double>();
		cms_scores_mean = new TreeMap<SNP, Double>();
		
		this.log = log;
	}
	
	public void runCmsAnalysis() throws FileParsingException {
		
		log.addLine("\nRunning Analysis");
		log.addLine("Importing Sim Data");
		
		System.out.println("\nRunning Analysis");
		System.out.println("Importing Sim Data");
		
		SimulationParser sp = new SimulationParser(log, sim_dir);
		SimDist[] neutral_sim = sp.getNeutralSimulations();
		SimDist[] select_sim = sp.getSelectedSimulations();
		
		log.addLine("Calculating Product and Mean of Scores");
		System.out.println("Calculating Product and Mean of Scores");
		for(int i = 0; i < all_ws.size(); i++) {
			WindowStats cur_ws = all_ws.get(i);
			
			System.out.println("Analyzing window with positions:\t" + cur_ws.getStPos() + " to " + cur_ws.getEndPos());//TODO: remove this
			calcCmsScores(cur_ws, neutral_sim, select_sim);
		}
		
		System.out.println("Analysis Finished");
	}
	
	public void writeStats() throws FileParsingException {
		
		System.out.println("Writing Data");
		
		try {
			out_dir = new File(out_dir.getAbsolutePath() + File.separator + "final_out");
			if(!out_dir.exists())
				out_dir.mkdir();
			
			File out_file = new File(out_dir.getAbsoluteFile() + File.separator 
					+ "analyzed_windows.tsv");
			int num = 1;
			while(out_file.exists()) {
				out_file = new File(out_dir.getAbsoluteFile() + File.separator 
						+ "analyzed_windows" + num + ".tsv");
				num++;
			}
			out_file.createNewFile();
			
			PrintWriter pw = new PrintWriter(out_file);
			pw.write(this.toString());
			pw.close();
			
			
		} catch(IOException e) {
			String msg = "Error: There was a problem with printing to the output file in " + out_dir.getAbsolutePath();
			throw new FileParsingException(log, msg);
		}
	}
	
	/*
	 * Simulations and scores need to be stored into arrays where:
	 * 		[0] = iHS data
	 * 		[1] = iHH data
	 * 		[2] = Fst data
	 * 		[3] = dDAF data
	 * 		[4] = XPEHH data
	 */
	private void calcCmsScores(WindowStats ws, SimDist[] neut_sim, SimDist[] sel_sim) 
			throws FileParsingException {
		
		TreeMap<SNP, Double> win_scores_prod = new TreeMap<SNP, Double>();
		TreeMap<SNP, Double> win_scores_mean = new TreeMap<SNP, Double>();
		
		double prior_prob = 1 / (double) ws.getTotNumSNPs();
		
		List<SNP> all_snps = ws.getAllSNPs();
		Double[] tst_scores = new Double[NUM_TESTS];
		
		for(int i = 0; i < all_snps.size(); i++) {
			
			SNP cur_snp = all_snps.get(i);
			System.out.println("\t" + cur_snp);
			
			tst_scores[0] = ws.getIhsScore(cur_snp);
			tst_scores[1] = ws.getIhhScore(cur_snp);
			tst_scores[2] = ws.getFstScore(cur_snp);
			tst_scores[3] = ws.getDDafScore(cur_snp);
			tst_scores[4] = ws.getXpehhScore(cur_snp);
			
			Double prod_cms = calcProductCMS(cur_snp, tst_scores, neut_sim, sel_sim, prior_prob);
			if(!prod_cms.equals(Double.NaN))
				win_scores_prod.put(cur_snp, prod_cms);
			
			Double mean_cms = calcMeanCMS(cur_snp, tst_scores, neut_sim, sel_sim, prior_prob);
			if(!mean_cms.equals(Double.NaN))
				win_scores_mean.put(cur_snp, mean_cms);
		}
		
		win_scores_prod = normalizeData(win_scores_prod);
		for(SNP key : win_scores_prod.keySet())
			cms_scores_prod.put(key, win_scores_prod.get(key));

		win_scores_mean = normalizeData(win_scores_mean);
		for(SNP key : win_scores_mean.keySet())
			cms_scores_mean.put(key, win_scores_mean.get(key));
		
		//TODO: analyze level of selection in entire window
	}
	
	private Double calcMeanCMS(SNP cur_snp, 
								Double[] tst_scores, 
								SimDist[] neut_sim, 
								SimDist[] sel_sim, 
								double prior_prob) throws FileParsingException {
		
		Double[] score_probs = new Double[NUM_TESTS];
		
		for(int j = 0; j < NUM_TESTS; j++) {
			
			if(!tst_scores[j].equals(Double.NaN)) {
			
				boolean two_sided = false;
				if(j == 0 || j== 1 || j == 3)
					two_sided = true;
				
				Double neut_prob = neut_sim[j].getProb(tst_scores[j], two_sided);
				Double sel_prob = sel_sim[j].getProb(tst_scores[j], two_sided);
				
				double cms_nom = sel_prob * prior_prob;
				double cms_denom = ((sel_prob*prior_prob) + (neut_prob*(1-prior_prob)));
				 
				score_probs[j] = cms_nom / cms_denom;
			}
		}
		
		//Mean of scores composition method
		Double final_score_mean = meanOfScores(score_probs);
		
		return final_score_mean;
	}
	
	private Double calcProductCMS(SNP cur_snp, 
									Double[] tst_scores, 
									SimDist[] neut_sim, 
									SimDist[] sel_sim, 
									double prior_prob) throws FileParsingException {
		
		Double[] score_probs = new Double[NUM_TESTS];
		
		boolean complete_data = true;
		for(int j = 0; j < NUM_TESTS; j++) {
			
			if(tst_scores[j].equals(Double.NaN) 
					|| getWS(cur_snp).getDafScore(cur_snp) < DAF_CUTOFF) {//Ilya told us to do the DAF_CUTOFF
				complete_data = false;
				break;
			}
			else {
			
				boolean two_sided = false;
				if(j == 0 || j== 1 || j == 3)
					two_sided = true;
				
				Double neut_prob = neut_sim[j].getProb(tst_scores[j], two_sided);
				Double sel_prob = sel_sim[j].getProb(tst_scores[j], two_sided);
				
				/*
				 * To do the CMSgw 
				 * 		I create SimulationParsers from different files
				 * 		Run a similar for loop structure
				 * 		instead of the below steps I do: double cms_bf (bayes factor) = sel_prob / neut_prob
				 * 		ask the question: Do I do this in parallel with the CMSlocal analysis
				 */
				
				double cms_nom = sel_prob * prior_prob;
				double cms_denom = ((sel_prob*prior_prob) + (neut_prob*(1-prior_prob)));
				 
				score_probs[j] = cms_nom / cms_denom;
			}
		}
		
		//Product of scores composition
		Double final_score = Double.NaN;
		if(complete_data)
			final_score = productOfScores(score_probs);
		
		return final_score;
	}
	
	private WindowStats getWS(SNP s) {
		for(WindowStats ws : all_ws) {
			if(ws.containsSNP(s))
				return ws;
		}
		
		return null;
	}
	
	@SuppressWarnings("unchecked")
	private TreeMap<SNP, Double> normalizeData(TreeMap<SNP, Double> unstd_cms) {
		
		List<Double> all_values = new LinkedList<Double>();
		TreeMap<SNP, Double> std_cms = new TreeMap<SNP, Double>();
		
		Iterator<SNP> itr = unstd_cms.navigableKeySet().iterator();
		while(itr.hasNext()) 
			all_values.add(unstd_cms.get(itr.next()));
		
		all_values = HaplotypeTests.normalizeData(all_values);
		
		itr = unstd_cms.navigableKeySet().iterator();
		int indx = 0;
		while(itr.hasNext()) {
			std_cms.put(itr.next(), all_values.get(indx));
			indx++;
		}
		
		return std_cms;
	}
	
	private Double productOfScores(Double[] score_probs) {
		
		Double score = 1.0;
		
		for(int i = 0; i < NUM_TESTS; i++) {
			if(score_probs[i] != null) 
				score = score*score_probs[i];
		}
		
		return score;
	}
	
	private Double meanOfScores(Double[] score_probs) {
		
		int tot_tests = NUM_TESTS;
		Double score = 0.0;
		for(int i = 0; i < NUM_TESTS; i++) {
			if(score_probs[i] != null)
				score += score_probs[i];
			else
				tot_tests--;
		}
		
		return score / tot_tests;
	}

	@Override
	public String toString() {
		
		List<SNP> all_keys = new LinkedList<SNP>();
		
		Iterator<SNP> itr = cms_scores_mean.navigableKeySet().iterator();
		while(itr.hasNext()) 
			all_keys.add(itr.next());
		
		StringBuilder sb = new StringBuilder();
		WindowStats ws = new WindowStats();
		sb.append("snp_id\tposition\tiHS\tXPEHH\tiHH\tdDAF\tDAF\tFst\tCMS_product\tCMS_mean\n");
		for(int i = 0; i < all_keys.size(); i++)  {
			SNP key = all_keys.get(i);
			
			if(!ws.containsSNP(key))
				ws = getWS(key);
			
			Double prod_cms = cms_scores_prod.get(key);
			if(prod_cms == null)
				prod_cms = Double.NaN;
			Double mean_cms = cms_scores_mean.get(key);
			if(mean_cms == null)
				mean_cms = Double.NaN;
			
			sb.append(key.getSnpID() + "\t" + key.getPosition() + "\t"
					+ ws.getIhsScore(key) + "\t" + ws.getXpehhScore(key) + "\t"
					+ ws.getIhhScore(key) + "\t" + ws.getDDafScore(key) + "\t"
					+ ws.getDafScore(key) + "\t" + ws.getFstScore(key) + "\t"
					+ prod_cms + "\t" + mean_cms + "\n");
		}
		
		return sb.toString();
	}

}
