from matplotlib import pyplot as plt
time_list = []
ecs_anomaly_list = []
true_anomaly_list = []
mean_anomaly_list = []

def read_file(name:str, output: list):
    with open(name,"r") as f:
        for _ in range(557):
            line_list = f.readline().split()
            time_list.append(int(line_list[0]))
            output.append(float(line_list[1]))

def ecs_graph():

    read_file("ecsAnomaly.txt", ecs_anomaly_list)
    plt.plot(time_list,ecs_anomaly_list)
    plt.xlabel("Время")
    plt.ylabel("Эксцентрическая аномалия")
    plt.show()

def mean_graph():

    read_file("orbit_data.txt", mean_anomaly_list)
    plt.plot(time_list,mean_anomaly_list)
    plt.xlabel("Время")
    plt.ylabel("Средняя аномалия")
    plt.show()

def true_graph():

    read_file("trueAnomaly.txt", true_anomaly_list)
    plt.plot(time_list,true_anomaly_list)
    plt.xlabel("Время")
    plt.ylabel("Истинная аномалия")
    plt.show()

true_graph()