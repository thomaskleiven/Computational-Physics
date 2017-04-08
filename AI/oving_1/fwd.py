import numpy as np
import heapq

transition_probability = np.matrix([[0.7,0.3], [0.3, 0.7]])
observationMatrix_true = np.matrix([[0.9, 0], [0,0.2]])
observationMatrix_false= np.matrix([[0.1, 0], [0, 0.8]])

class fwd_bwd:
    def fwd(ev, pri):
        fv = np.array([None]*(len(ev)+1))          #A vector of forward messages for teps 0,...,t
        fv[0] = pri
        alphas = np.zeros(len(ev))
        print("Running forward: ")
        print("Message 0: ", fv[0])
        for i in range(1,len(ev)+1):
            if(ev[i-1] == True):
                fwd = np.dot(observationMatrix_true*np.transpose(transition_probability), fv[i-1])
            else:
                fwd = np.dot(observationMatrix_false*np.transpose(transition_probability), fv[i-1])
            fv[i] = fwd_bwd.normalize(np.ravel(fwd.sum(axis=0)))
            print("Message %d: %s"%(i,fv[i]))
        return fv

    def normalize(vec):
        return vec/vec.sum()

    def bwd(ev, fv):
        b = np.array([1,1])                         #A representation of the backward messages, initially all 1s
        sv = np.array([None]*(len(ev)+1))           #A vector of smoothed estimates for steps 1,...,t
        print(" ")
        print("Running backwards: ")
        print("Message 0: ", b)

        t = len(ev)+1
        for i in range(t-1, -1, -1):
            print(fv[i])
            print(b)

            sv[i] = fwd_bwd.normalize(fv[i] * b)
            print(sv[i])
            if(ev[i-1] == True):
                b = np.ravel((np.dot(transition_probability*observationMatrix_true,b)).sum(axis=0))
            else:
                b = np.ravel((np.dot(transition_probability*observationMatrix_false,b)).sum(axis=0))
            print("Message %d: %s"%(i, b))
        return sv

    def viterbi(fv, ev):
        sequence = []
        temporaryGuess = np.array([None]*(len(fv)+1))
        temporaryGuess[0] = fv[0]
        b = np.array(0)
        for i in range(0,len(fv)):
            if(ev[i] == True):
                a = transition_probability*observationMatrix_true
                b = (a[:,:1])
                c = np.multiply(temporaryGuess[i],np.transpose(b))

                d = (np.ravel(a)[1::2])
                e = d.reshape(2,1)
                f = np.multiply(temporaryGuess[i], np.transpose(e))

                temporaryGuess[i+1] = max(np.ravel(c)), max(np.ravel(f))
            else:
                a = transition_probability*observationMatrix_false
                b = (a[:,:1])
                c = np.multiply(temporaryGuess[i],np.transpose(b))

                d = (np.ravel(a)[1::2])
                e = d.reshape(2,1)
                f = np.multiply(temporaryGuess[i], np.transpose(e))

                temporaryGuess[i+1] = max(np.ravel(c)), max(np.ravel(f))
        return temporaryGuess


def main():
    evidence = np.array([True, True, False, True, True])
    prior = np.array([0.5, 0.5])
    result = fwd_bwd.fwd(evidence, prior)
    back = fwd_bwd.bwd(evidence, result)
    print("--------------------------")
    print("Results: ")
    for i in range(0,len(back)):                     #Print results
        print("%d: %s"%(i, back[i]))
    print("--------------------------")
    print("Results Viterbi: ")
    vit = fwd_bwd.viterbi(result[1:], evidence)
    for i in range(1,len(vit)):
        if(vit[i][0] > vit[i][1]):
            print("True")
        else:
            print("False")

main()
