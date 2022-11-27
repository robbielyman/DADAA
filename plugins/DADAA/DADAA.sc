TADAA : MultiOutUGen {
	*ar { |input, gain|
		^this.multiNew('audio', input, gain);
	}
  init { arg ... theInputs;
    inputs = theInputs;
    ^this.initOutputs(2, rate);
  }
	checkInputs {
		^this.checkValidInputs;
	}
}

CADAA : MultiOutUGen {
	*ar { |input, gain|
		^this.multiNew('audio', input, gain);
	}
  init { arg ... theInputs;
    inputs = theInputs;
    ^this.initOutputs(2, rate);
  }
	checkInputs {
		^this.checkValidInputs;
	}
}
