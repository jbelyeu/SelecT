package stats;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.PrintWriter;
import java.util.List;

import calc.*;
import tools.*;

public class StatsCalc {

	/** Hyper-Parallelized Composite of Multiple Signals (CMS) Java implementation GWAS and Local study version 1.0
	 * This program calculated the stats that are used in CMS analysis
	 * @author Hayden Smith
	 * 
	 * @param args0		out dir
	 * @param args1		chr
	 * @param args2		window number
	 */
	public static void main(String[] args) {
		
		File out_dir = new File(args[0]);
		if(!out_dir.exists() && !out_dir.isDirectory()) {
			
			Log log = new Log(Log.type.stat);
			log.addLine("ERROR: Parameter for out file directory invalid");
			log.addLine("\t*Go to api for more information");
			log.addLine("\t*You will need to redo this entire step--all new data is invalid");
			log.close();
			
			System.exit(0);
		}
		
		int chr = -1;
		int win_num = -1;
		try {
			
			chr = Integer.parseInt(args[1]);
			win_num = Integer.parseInt(args[2]);
			
		} catch(NumberFormatException e) {
			
			Log log = new Log(Log.type.stat);
			log.addLine("ERROR: Parameters for chr and/or window number invalid");
			log.addLine("\t*Go to api for more information");
			log.addLine("\t*You will need to redo this entire step--all new data is invalid");
			log.close();
			
			System.exit(0);
		}
		
		Log log = new Log(Log.type.stat, out_dir.getName());
		try {
			System.out.println("Running Window: " + win_num);
			System.out.println("Chromosome: " + chr);
			System.out.println("Output Dir:" + out_dir);
			
			StatsCalc sc = new StatsCalc(out_dir, chr, win_num, log);
			sc.runStats();
			
			System.out.println("Window complete!");
			log.addLine(win_num + "\tSuccessfulRun\twindow completed without any errors");
			log.close();
			
		 } catch(OutOfMemoryError e) {
			 
			 log.addLine(win_num + "\tNoMemoryError\tinsufficient memory for running this window");
			 log.close();
			 
			 System.exit(0);
		 }
	}
	
	private static int WAIT_TIME = 50;
	
	private int chr;
	private int win_num;
	
	private File out_dir;
	private File envi_dir;
	private File win_dir;
	
	private Window tp_win;
	
	private iHS i;
	private iHH h;
	private XPEHH x;
	private dDAF d;
	private Fst f;
	//private TajD t;
	//private NewStat new;
	
	private Log log;
	
	public StatsCalc(File out_dir, int chr, int win_num, Log log) {
		
		tp_win = null;
		
		this.chr = chr;
		this.win_num = win_num;
		this.out_dir = out_dir;
		this.log = log;
	}
	
	public void runStats() {
		
		setupFiles();
		
		createCalculators();
		
		doCalculations();
		
		writeOutput();
	}
	
	private void writeOutput() {
		
		try {
			
			File win_stats = new File(out_dir.getAbsoluteFile() + File.separator 
					+ "win" + win_num + "_" + "chr" + chr + "_s" 
					+ tp_win.getStPos() + "-e" + tp_win.getEndPos() + ".tsv");
			win_stats.createNewFile();
			
			WindowStats ws = new WindowStats(tp_win.getStPos(), tp_win.getEndPos());
			
			ws.setIHS(i.getStats(), i.getSNPs());
			ws.setIHH(h.getStats(), h.getSNPs());
			ws.setXPEHH(x.getStats(), x.getSNPs());
			ws.setDDAF(d.getStats(), d.getSNPs());
			ws.setDAF(d.getDafStats(), d.getSNPs());
			ws.setFst(f.getStats(), f.getSNPs());
			//ws.setTAJD(t.getStats(), t.getSNPs());
			//ws.setNEW(new.getStats(), new.getSNPs());
			
			PrintWriter pw = new PrintWriter(win_stats);
			pw.print("snp_id\tposition\tiHS\tXPEHH\tiHH\tdDAF\tDAF\tFst\n");//\tTajD\tNew
			pw.print(ws);
			pw.close();
			
		} catch (IOException e) {
			log.addLine(win_num + "\tWriteFileError\tCould not create/write the output file for this window");
			System.exit(0);
		}
	}
	
	private void doCalculations() {
			
		Object lock = new Object();
		
		try {
		
			StatsThread i_thrd = new StatsThread(i, lock);
			Thread.sleep(WAIT_TIME);
			i_thrd.start();
			
			StatsThread h_thrd = new StatsThread(h, lock);
			Thread.sleep(WAIT_TIME);
			h_thrd.start();
			
			StatsThread x_thrd = new StatsThread(x, lock);
			Thread.sleep(WAIT_TIME);
			x_thrd.start();
			
			StatsThread d_thrd = new StatsThread(d, lock);
			Thread.sleep(WAIT_TIME);
			d_thrd.start();
			
			StatsThread f_thrd = new StatsThread(f, lock); 
			Thread.sleep(WAIT_TIME);
			f_thrd.start();
			
			//StatsThread t_thrd = new StatsThread(t, lock);
			//Thread.sleep(WAIT_TIME);
			//t_thrd.start();
			
			//StatsThread new_thrd = new StatsThread(new, lock);
			//Thread.sleep(WAIT_TIME);
			//new_thrd.start();
			
			synchronize(i_thrd, h_thrd, x_thrd, d_thrd, f_thrd);//add t_thrd
			
			i = (iHS) i_thrd.getTest();
			h = (iHH) h_thrd.getTest();
			x = (XPEHH) x_thrd.getTest();
			d = (dDAF) d_thrd.getTest();
			f = (Fst) f_thrd.getTest();
			//t = (TajD) t_thrd.getTest();
			//new = (NewStat) new_thrd.getTest();
		
		} catch (InterruptedException e) { 
			log.addLine(win_num + "\tThreadingError\tThreading for this process did not complete");
			System.exit(0);
		}
	}
	
	private void synchronize(StatsThread i_thrd, 
								StatsThread h_thrd, 
								StatsThread x_thrd, 
								StatsThread d_thrd, 
								StatsThread f_thrd) {

		for(;;) {
			if(i_thrd.isFinished()
					&& h_thrd.isFinished()
					&& x_thrd.isFinished()
					&& d_thrd.isFinished()
					&& f_thrd.isFinished())
				break;
			else
				continue;
		}
	}
	
	@SuppressWarnings("unchecked")
	private void createCalculators() {
		
		try {
			String path = "";
			
			path = getEnviFileName("target_pop_wins.bin");
			List<Window> tp_wins = (List<Window>) getObject(path);
			
			path = getEnviFileName("target_pop_indv.bin");
			Individual[] tp_indv = (Individual[]) getObject(path);
			
			path = getEnviFileName("targetXcross_wins.bin");
			List<Window> txin_wins = (List<Window>) getObject(path);
			
			path = getEnviFileName("targetXcross_indv.bin");
			Individual[] tp_inx_indv = (Individual[]) getObject(path);
			
			path = getEnviFileName("crossXtarget_indv.bin");
			Individual[] xp_int_indv = (Individual[]) getObject(path);
			
			path = getEnviFileName("crossXout_wins.bin");
			List<Window> xoin_wins = (List<Window>) getObject(path);
			
			path = getEnviFileName("crossXout_indv.bin");
			Individual[] xp_ino_indv = (Individual[]) getObject(path);
			
			path = getEnviFileName("outXcross_indv.bin");
			Individual[] op_inx_indv = (Individual[]) getObject(path);
			
			path = getEnviFileName("anc_types.bin");
			List<Window> anc_types = (List<Window>) getObject(path);
			
			path = getEnviFileName("genetic_map.bin");
			GeneticMap gm = (GeneticMap) getObject(path);
			
			path = getTargetWindowFileName();
			tp_win = (Window) getObject(path);
			
			path = getTargetXCrossWindowFileName();
			Window txin_win = (Window) getObject(path);
			
			i = new iHS(tp_win, tp_indv, anc_types, tp_wins, gm);
			h = new iHH(tp_win, tp_indv, anc_types, tp_wins, gm);
			x = new XPEHH(txin_win, txin_wins, tp_inx_indv, xp_int_indv, gm);
			d = new dDAF(tp_win, tp_indv, xoin_wins, xp_ino_indv, op_inx_indv, anc_types);
			f = new Fst(txin_win, tp_inx_indv, xp_int_indv, op_inx_indv);
			//t = new TajD(args0, args1, ..., argsx);
			//new = new NewStat(args0, args1, ..., argsx);
			
		} catch (IOException e) {
			log.addLine(win_num + "\tReadFileError\tCould not find the correct file for proper loading of envi; check chr num and api");
			System.exit(0);
		} catch (ClassNotFoundException e) {
			log.addLine(win_num + "\tClassNotFoundError\tCould create the correct object instance while loading of envi");
			System.exit(0);
		} catch (ClassCastException e) {
			log.addLine(win_num + "\tCastingError\tObject casting invalid while loading envi");
			System.exit(0);
		}
	}
	
	@SuppressWarnings("resource")
	private Object getObject(String path) throws IOException, ClassNotFoundException {
		
		File file = new File(path);
		if(!file.exists()) {
			throw new IOException();
		}
		
		ObjectInputStream ois = new ObjectInputStream(new BufferedInputStream(new FileInputStream(file)));
		return ois.readObject();
	}
	
	private String getTargetXCrossWindowFileName() throws IOException {
		
		String[] all_files = win_dir.list();
		for(int i = 0; i < all_files.length; i++) {
			
			String file_name = all_files[i];
			if(file_name.contains("chr" + chr)
					&& file_name.contains("win" + win_num + "_")
					&& file_name.contains("x"))
				return win_dir.getAbsolutePath() + File.separator + file_name;
		}
		
		throw new IOException();
	}
	
	private String getTargetWindowFileName() throws IOException {
		
		String[] all_files = win_dir.list();
		for(int i = 0; i < all_files.length; i++) {
			
			String file_name = all_files[i];
			if(file_name.contains("chr" + chr)
					&& file_name.contains("win" + win_num)
					&& !file_name.contains("x")) 
				return win_dir.getAbsolutePath() + File.separator + file_name;
		}
		
		throw new IOException();
	}
	
	private String getEnviFileName(String name) {
		return envi_dir.getAbsolutePath() + File.separator + name;
	}
	
	private void setupFiles() {
		
		envi_dir = new File(out_dir.getAbsoluteFile() + File.separator + "envi_files" + File.separator + "envi_var");
		if(!envi_dir.exists() && !envi_dir.isDirectory()) {
			log.addLine(win_num + "\tEnviDirError\tCould not find the evironment directory; check api for path names");
			System.exit(0);
		}
		
		win_dir = new File(out_dir.getAbsoluteFile() + File.separator + "envi_files" + File.separator + "all_wins");
		if(!win_dir.exists() && !win_dir.isDirectory()) {
			log.addLine(win_num + "\tWindowDirError\tCould not find the windows directory; check api for path names");
			System.exit(0);
		}
		
		out_dir = new File(out_dir.getAbsolutePath() + File.separator + "stats_files");
		if(!out_dir.exists())
			out_dir.mkdirs();
	}
}

class StatsThread extends Thread {

	private final Object lock;
	
	private Thread thrd;
	private HaplotypeTests tst;
	
	volatile private boolean finished; //volatile says that some other thread could change this value
	
	StatsThread(HaplotypeTests tst, Object lock) {
		this.tst = tst;
		this.lock = lock;
		
		finished = false;
		
		thrd = new Thread(this);
//		thrd.start();
	}
	
	public void start() {
		thrd.start();
	}
	
	@Override
	public void run() {
		
		tst.runStat();	
		
		synchronized(lock) {
			finished = true;
		}
		
		thrd.interrupt();
	}
	
	public HaplotypeTests getTest() {
		return tst;
	}
	
	public boolean isFinished() {
		
		synchronized(lock) {
			return finished;
		}
	}
	
}
