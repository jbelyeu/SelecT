package analysis;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.LinkedList;
import java.util.List;
import java.util.Scanner;

import errors.FileParsingException;
import tools.*;

public class WindowParser {
	
	private File win_file;
	
	private WindowStats ws;
	
	private Log log;
	
	public WindowParser(Log log, File win_file, int st_pos, int end_pos) {
		
		this.log = log;
		this.win_file = win_file;
		
		ws = new WindowStats(st_pos, end_pos);
	}
	
	public WindowStats parseWindow() throws FileParsingException {
		
		List<SNP> ihs_snps = new LinkedList<SNP>();
		List<SNP> xpehh_snps = new LinkedList<SNP>();
		List<SNP> ihh_snps = new LinkedList<SNP>();
		List<SNP> ddaf_snps = new LinkedList<SNP>();
		List<SNP> Xdaf_snps = new LinkedList<SNP>();
		List<SNP> fst_snps = new LinkedList<SNP>();
		//List<SNP> tajd_snps = new LinkedList<SNP>();
		//List<SNP> new_snps = new LinkedList<SNP>();
		
		List<Double> ihs_stats = new LinkedList<Double>();
		List<Double> xpehh_stats = new LinkedList<Double>();
		List<Double> ihh_stats = new LinkedList<Double>();
		List<Double> ddaf_stats = new LinkedList<Double>();
		List<Double> Xdaf_stats = new LinkedList<Double>();
		List<Double> fst_stats = new LinkedList<Double>();
		//List<Double> fajd_stats = new LinkedList<Double>();
		//List<Double> new_stats = new LinkedList<Double>();
		
		try {
			
			Scanner scan = new Scanner(win_file);
			scan.nextLine();//skip header line
			
			while(scan.hasNext()) {
				
				String[] line = scan.nextLine().split("\\s+");
				
				if(line.length != 8) {
					String msg = "Error: Window file " + win_file.getName() + " has irregular number of scores";
					throw new FileParsingException(log, msg);
				}
				
				SNP s = new SNP(Integer.parseInt(line[1]), line[0]);
				Double i_dbl = Double.parseDouble(line[2]);
				Double x_dbl = Double.parseDouble(line[3]);
				Double h_dbl = Double.parseDouble(line[4]);
				Double dd_dbl = Double.parseDouble(line[5]);
				Double d_dbl = Double.parseDouble(line[6]);
				Double f_dbl = Double.parseDouble(line[7]);
				//Double t_dbl = Double.parseDouble(line[8]);
				//Double new_dbl = Double.parseDouble(line[X]);
				
				if(i_dbl != Double.NaN) {
					ihs_snps.add(s);
					ihs_stats.add(i_dbl);
				}
				if(x_dbl != Double.NaN) {
					xpehh_snps.add(s);
					xpehh_stats.add(x_dbl);
				}
				if(h_dbl != Double.NaN) {
					ihh_snps.add(s);
					ihh_stats.add(h_dbl);
				}
				if(dd_dbl != Double.NaN) {
					ddaf_snps.add(s);
					ddaf_stats.add(dd_dbl);
				}
				if(d_dbl != Double.NaN) {
					Xdaf_snps.add(s);
					Xdaf_stats.add(d_dbl);
				}
				if(f_dbl != Double.NaN) {
					fst_snps.add(s);
					fst_stats.add(f_dbl);
				}
				//if(t_dbl != Double.NaN) {
				//	tajd_snps.add(s);
				//	tajd_stats.add(t_dbl);
				//}
				//if(new_dbl != Double.NaN) {
				//	new_snps.add(s);
				//	new_stats.add(new_dbl);
				//}
			}
			
		} catch (FileNotFoundException e) {
			String msg = "Error: Could not open data file";
			throw new FileParsingException(log, msg);
		} catch (NumberFormatException e) {
			String msg = "Error: Irregular data types in window windo file " + win_file.getName();
			throw new FileParsingException(log, msg);
		}
		
		ws.setIHS(ihs_stats, ihs_snps);
		ws.setIHH(ihh_stats, ihh_snps);
		ws.setXPEHH(xpehh_stats, xpehh_snps);
		ws.setDDAF(ddaf_stats, ddaf_snps);
		ws.setDAF(Xdaf_stats, Xdaf_snps);
		ws.setFst(fst_stats, fst_snps);
		//ws.setTAJD(tajd_stats, tajd_snps);
		//ws.setNEW(new_stats, new_snps);
		
		return ws;
	}

}
