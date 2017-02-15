import numpy as np

transition_probability = np.matrix([[0.7,0.3], [0.3, 0.7]])
observationMatrix_true = np.matrix([[0.9, 0], [0,0.2]])
observationMatrix_false= np.matrix([[0.1, 0], [0, 0.8]])

initial_condition = np.array([0.5, 0.5])

class fwd_bwd:
    def fwd(ev, pri):
        fv = np.array(len(ev)) #A vector of forward messages for teps 0,...,t
        fv = pri
        alphas = np.zeros(len(ev))
        for i in range(0,len(ev)):
            total_probability = fv[0]*np.transpose(transition_probability)[0] + fv[1]*np.transpose(transition_probability)[1]
            if(ev[i] == True):
                bayes1 = observationMatrix_true[0]*np.transpose(total_probability)
                bayes2 = observationMatrix_true[1]*np.transpose(total_probability)
            else:
                bayes1 = observationMatrix_false[0]*np.transpose(total_probability)
                bayes2 = observationMatrix_false[1]*np.transpose(total_probability)
            alpha = 1/(bayes1 + bayes2)
            alphas[i] = alpha
            fv[0] = bayes1*alpha
            fv[1] = bayes2*alpha
        return fv

    def bwd(ev, fv):
        sv = np.array(len(ev)) #A smoothed estimate for steps 1,...,t
        b = np.ones(len(ev)) #A vector of backward messages
        for i in range(0, len(ev)):
            total_probability = fv[i]*b[i]
            b = transition_probability*observationMatrix_true
        print(observationMatrix_true)
        print(transition_probability)
        print(np.transpose(b))



def main():
    evidence = np.array([True, True])
    prior = np.array([0.5, 0.5])
    result = fwd_bwd.bwd(evidence, fwd_bwd.fwd(evidence, prior))

main()
