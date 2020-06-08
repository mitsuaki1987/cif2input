#!/usr/bin/python3
import combo
import numpy
import sys


def load_data(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
        result = []
        descriptor = []
        for line in lines:
            result.append(float(line.split()[0]))
            descriptor.append(numpy.array(line.split()[1:], dtype=numpy.float_))

    return numpy.array(descriptor), numpy.array(result)


class Simulator:
    def __init__(self, filename):
        _, self.result = load_data(filename)

    def __call__(self, action):
        return self.result[action]


def main():
    args = sys.argv
    if len(args) < 2:
        print("Usage:")
        print("$ mi_test.py result-desc-file [rand_seed] [max_action] [max_rand]")
        print("Default:")
        print("$ cif2input.py result-desc-file 1 200 10")
        exit(0)

    rand_seed = 1
    max_action = 200
    max_rand = 10
    result_desc = args[1]
    if len(args) > 2:
        rand_seed = int(args[2])
        if len(args) > 3:
            max_action = int(args[3])
            if len(args) > 4:
                max_rand = int(args[4])

    descriptor, exact = load_data(result_desc)
    simulator = Simulator(result_desc)
    descriptor = combo.misc.centering(descriptor)
    policy = combo.search.discrete.policy(test_X=descriptor)
    policy.set_seed(rand_seed)
    #
    f = open("history.dat", 'w')
    for i_action in range(max_action):
        if i_action < max_rand:
            action = policy.random_search(max_num_probes=1, num_search_each_probe=1, simulator=None)
        else:
            action = policy.bayes_search(max_num_probes=1, num_search_each_probe=1, simulator=None,
                                         score='TS', interval=0, num_rand_basis=0)
        result = simulator(action)
        if result < 1.0e-5:
            policy.delete_actions(action)
        else:
            policy.write(action, result)
        print(i_action, action[0], result[0],
              numpy.max(policy.history.fx[0:policy.history.total_num_search]),
              file=f, flush=True
              )
    f.close()
    #
    # Print regression and error
    #
    X_train = numpy.array(descriptor[policy.history.chosed_actions[0:policy.history.total_num_search], :])
    y_train = numpy.array(exact[policy.history.chosed_actions[0:policy.history.total_num_search]])

    cov = combo.gp.cov.gauss(X_train.shape[1], ard=False)
    mean = combo.gp.mean.const()
    lik = combo.gp.lik.gauss()
    gp = combo.gp.model(lik=lik, mean=mean, cov=cov)
    config = combo.misc.set_config()
    gp.fit(X_train, y_train, config)
    gp.prepare(X_train, y_train)
    fmean = gp.get_post_fmean(X_train, descriptor)
    fcov = gp.get_post_fcov(X_train, descriptor)

    with open("combo.dat", mode="w") as f:
        for i in range(len(exact)):
            print(exact[i], fmean[i], numpy.sqrt(fcov[i]), file=f)

    fmean = gp.get_post_fmean(X_train, X_train)
    fcov = gp.get_post_fcov(X_train, X_train)

    with open("combo2.dat", mode="w") as f:
        for i in range(len(y_train)):
            print(y_train[i], fmean[i], numpy.sqrt(fcov[i]), file=f)


main()
