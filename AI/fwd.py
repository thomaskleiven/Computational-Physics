import numpy as np

transition_probability = np.matrix([[0.7,0.3], [0.3, 0.7]])
observationMatrix_true = np.matrix([[0.9, 0], [0,0.2]])
observationMatrix_false= np.matrix([[0.1, 0], [0, 0.8]])

initial_condition = np.array([0.5, 0.5])

class fwd_bwd:
    def fwd(ev, pri):
        fv = np.array([None]*len(ev)) #A vector of forward messages for teps 0,...,t
        fv[0] = pri
        alphas = np.zeros(len(ev))
        for i in range(1,len(ev)):
            if(ev[i] == True):
                fwd = np.dot(observationMatrix_true*np.transpose(transition_probability), fv[i-1])
            else:
                fwd = np.dot(observationMatrix_false*np.transpose(transition_probability), fv[i-1])
            fv[i] = fwd_bwd.normalize(np.ravel(fwd.sum(axis=0)))

    def normalize(vec):
        return vec/vec.sum()    


def main():
    evidence = np.array([True, True, False, True, True])
    prior = np.array([0.5, 0.5])
    fwd_bwd.fwd(evidence, prior)
main()
