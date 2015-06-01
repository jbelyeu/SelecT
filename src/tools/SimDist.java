package tools;

import java.util.ArrayList;
import java.util.List;

import errors.FileParsingException;

public class SimDist {
	
	private final int BIN_NUM = 60;
	
	private int up_bndry;
	private int low_bndry;
	private double total_prob;
	
	private List<Double> sim_vals;
	
	public Log log;
	
	public SimDist(Log log, int low_bndry, int up_bndry) {
		
		this.log = log;
		
		this.up_bndry = up_bndry;
		this.low_bndry = low_bndry;
		this.total_prob = 0.0;
		
		sim_vals = new ArrayList<Double>();
	}

	public void addSimValue(double val) {
		
		sim_vals.add(val);
		total_prob += val;
	}
	
	public Double get1SidedProb(Double score) throws FileParsingException {
		
		if(sim_vals.size() != BIN_NUM) {
			String msg = "Bin numbers don't coincide, error in reading simulated data";
			throw new FileParsingException(log, msg);
		}
		
		int s_indx = getScoreIndex(score, up_bndry, low_bndry);
		
		return calcProbAtBin(s_indx);
	}
	
	public Double get2SidedProb(Double score) throws FileParsingException {
		
		if(sim_vals.size() != BIN_NUM) {
			String msg = "Bin numbers don't coincide, error in reading simulated data";
			throw new FileParsingException(log, msg);
		}
		
		score = Math.abs(score);
		
		int up_indx = getScoreIndex(score, up_bndry, low_bndry);
		int low_indx = getScoreIndex((-1 * score), up_bndry, low_bndry);
		
		return calcTwoSidedProbAtBin(up_indx, low_indx);
	}
	
	/**
	 * Handles probability calculations for simulations. Uses inputed boolean to
	 * decided which kind of probability to run
	 * 
	 * @param score			Probability of this score given the SimDist score 
	 * 						distribution
	 * @param two_sided		If true this method calculates two-sided probability;
	 * 						if false it calls one-sided getProb function
	 * @return
	 * @throws FileParsingException
	 */
	public Double getProb(Double score, boolean two_sided) throws FileParsingException {
		
		if(two_sided)
			return get2SidedProb(score);
		else
			return get1SidedProb(score);
	}
	
	public List<Double> getSimVals() {
		return sim_vals;
	}
	
	public Double calcTotalProb() {
		
		Double prob = 0.0;
		
		for(int i = 0; i < sim_vals.size(); i++)
			prob += sim_vals.get(i);
		
		return prob;
	}
	
	public double getTotalProb() {
		return total_prob;
	}
	
	private Double calcProbAtBin(int indx) {
		
		Double prob = 0.0;
		
		for(int i = 0; i <= indx; i++)
			prob += sim_vals.get(i);
		
		return prob;
	}
	
	private Double calcTwoSidedProbAtBin(int up_indx, int low_indx)  {
		
		Double prob = 0.0;
		
		for(int i = low_indx; i <= up_indx; i++)
			prob += sim_vals.get(i);
		
		return prob;
	}
	
	private int getScoreIndex(Double score, int up, int dwn) {
		
		double rng = (double) up - (double) dwn;
		double bin_size = rng / (double) BIN_NUM;
		
		for(int i = 0; i < BIN_NUM; i++) {
			if((dwn + bin_size*i) >= score)
				return i;
		}
		
		return BIN_NUM-1;
	}
	
}
