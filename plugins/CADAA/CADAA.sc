CADAA : MultiOutUGen {
	*ar { |input|
		^this.multiNew('audio', input);
	}
  init { arg ... theInputs;
    inputs = theInputs;
    ^this.initOutputs(2, rate);
  }
	checkInputs {
		^this.checkValidInputs;
	}
}
