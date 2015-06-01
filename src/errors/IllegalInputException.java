package errors;

import tools.Log;

@SuppressWarnings("serial")
public class IllegalInputException extends Exception {
	
	public IllegalInputException() {}
	
	public IllegalInputException(Log log, String message) {
		log.addLine("\n");
		log.addLine(message);
		log.addLine("\t*Make sure all arguments are correct" );
		log.addLine("\t*Rerun with -h option as first argument in argument list for help");
		log.addLine("\t*Go to api for parameter descriptions");
	}

}
