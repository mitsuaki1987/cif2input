#!/usr/bin/python3
import combo
import numpy


def load_data():
    with open("mi.dat", "r") as f:
        lines = f.readlines()
        result = []
        descriptor = []
        for line in lines:
            result.append(float(line.split()[0]))
            descriptor.append(numpy.array(line.split()[1:], dtype=numpy.float_))

    return numpy.array(descriptor), numpy.array(result)


class Simulator:
    def __init__(self):
        _, self.result = load_data()

    def __call__(self, action):
        return self.result[action]


def main():
    descriptor, _ = load_data()
    simulator = Simulator()
    descriptor = combo.misc.centering(descriptor)
    policy = combo.search.discrete.policy(test_X=descriptor)
    policy.set_seed(3)
    #
    f = open("history.dat", 'w')
    max_action = 50
    for i_action in range(max_action):
        if i_action < 10:
            action = policy.random_search(max_num_probes=1, num_search_each_probe=1, simulator=None)
        else:
            action = policy.bayes_search(max_num_probes=1, num_search_each_probe=1, simulator=None,
                                         score='EI', interval=0, num_rand_basis=0)
        result = simulator(action)
        if result < 1.0e-5:
            policy.delete_actions(action)
        else:
            policy.write(action, result)
        print(i_action, action[0], result[0],
              # numpy.max(policy.history.fx[0:policy.history.total_num_search]),
              file=f
              )
    f.close()
    print(numpy.sort(policy.history.chosed_actions[0:policy.history.total_num_search]))


main()
