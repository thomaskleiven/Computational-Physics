import numpy as np

def main():
  # Initialize the variables for the umbrella domain given in chapter 15 of AIMA.
  O = np.array([[0.9, 0.0],[0.0, 0.2]])
  T = np.array([[0.7, 0.3],[0.3, 0.7]])
  prior = np.array([0.5, 0.5])
  evidence = [True, True, False, True, True]
  res = forward_backward(evidence, prior, O, T)
  print ("---- Results ---")
  for i in range(len(res)):
    print ("%i: %s" % (i, res[i]))

def forward_backward(ev, prior, O, T):
  """
    Forward-backward algorithm defined in figure 15.4 in AIMA.
    Computes probabilities of states given a sequence of obeservations (evidence)
    Input:
      - ev: this is the evidence and equals our observations
      - prior: the initial probabilities for the different states
      - T: transition matrix, calculates the ...
      - O:  Sensor matrix - probability for evidence variables given X
    Output:
      Given evidence values it then calculates and returns the probabilties
      for these evidences to be 'facts'. Returned in a array with distrubutions
      such as [[0.5,0.5],[0.883,0.117]]
  """

  def forward(fv, ev, O, T):
    """
     Input:
        fv: forward messages
        ev: evidence variable for i, e.g. True
        O:  Sensor matrix - probability for evidence variables given X
        T:  transition matrix

      Output:
        A vector of normalized probabilities, e.g. [0.818, 0.182]
    """
    if not ev:
      # Subtract from the identification matrix
      O = np.eye(len(O)) - O
    # Do calculations as in Equation 15.12 in AIMA
    r = np.dot(O, np.transpose(T))
    r = np.dot(r, fv)
    return normalize(r)

  def normalize(v):
    """
      Normalizes vector so that the sum/product equals to 1
      Input:
        v: A vector to be normalized
      Output:
        normalize([0.45, 0.10]) => [0.8182, 0.1818]
    """
    return v/v.sum()

  def backward(b, ev, O, T):
    """
      We want to use the results from forward in order to do a smoothing/better
      prediction of the probabilties. Therefore we complete the forward-backward
      algorithm with a backward step which takes the results from forward into
      account and do a matrix multuplication on this, the evidence and previous
      results.
      Input:
        b: Previous backward message, initially [1,1] (see forward_backward()).
        ev: The evidence state (this is what we are observing)
        O: Sensor matrix - probability for evidence variables given X
        T: Transition matrix
      Output:
        A new b (backward message), e.g. [0.690, 0.410]
    """
    if not ev:
      # Subtract from the identification matrix
      O = np.eye(len(O)) - O
    # Do calculations as in Equation 15.13 in AIMA
    r = np.dot(T, O)
    r = np.dot(r, b)
    return r

  # Initialize the forward messages array
  t = len(ev) + 1
  fv = np.array([None]*t)
  fv[0] = prior
  print ("---- Running forward ---")
  print ("Forward-message 0: %s" % fv[0])
  for i in range(1, t):
    # For each evidence variable we update the forward messages
    fv[i] = forward(fv[i - 1], ev[i-1], O, T)
    print ("Forward-message %i: %s" % (i, fv[i]))
  # Inititalize the smoothing and backwards array
  sv = np.array([None]*t)
  sv[0] = prior
  b = np.array([1,1])

  print ("---- Running backward ---")
  for j in range(t-1, -1, -1):
    # For each iteration we update the backward message variable
    # We then calculate the new probabilty, taken the old into account
    sv[j] = normalize( fv[j] * b )
    #print "Smoothing array for step %i: %s" % (j, sv[j])
    print ("Backward-array for step %i: %s " % (j, b))
    b = backward(b, ev[j-1], O, T)

  # Return the final result after smoothing
  return sv

if __name__ == '__main__':
  # Start the program
  main()
