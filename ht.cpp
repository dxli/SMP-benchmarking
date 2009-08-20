//test efficiency of multi-core/hyper-threading
// measure speed of each thread at various numbers of running processes
#include<iostream>
#include<fstream>
#include<vector>
#include<cstdlib>
#include<string>
#include<unistd.h>
#include<grace_np.h>
#include<sys/types.h>
#include<sys/wait.h>
#include<sys/time.h>
using namespace std;
int main(int argc,char *argv[])
{
        int np=4;
        if(argc>=2){
                np=atoi(argv[1]);
                if ( !( np >=2 && np <=32) ) np =4;
        }
        cout<<"testing up to "<<np<<" processes\n";
    vector<pid_t> spid;
    double speed0;
    vector<double> xi,yi;
	string fn("ht-benchmark.txt");
cout<<"writing to "<<fn<<endl;
	unlink(fn.c_str());
    for (int i=1;i<=np;i++) {
        struct timeval t_start,t_end;
        gettimeofday(&t_start,NULL);
        for (int j=1;j<=i;j++) {
            pid_t pid0=fork();
            if ( ! pid0 ) {
                execl("./b", "b", "21",NULL);
		exit(0);
            } else spid.push_back(pid0);
        }
        for (unsigned int j1=0;j1<spid.size();j1++) {
            int status;
            waitpid(-1,&status,WUNTRACED|WCONTINUED);
        }
	spid.clear();
        gettimeofday(&t_end,NULL);
        double t=(t_end.tv_sec - t_start.tv_sec) + 1e-6*(t_end.tv_usec - t_start.tv_usec);
        t_start=t_end;
        double speed= 1./0.114532/t;
        //if (i==1) speed0=speed;
	cout<<i<<' '<<speed<<endl;
        xi.push_back(double(i));
        yi.push_back(speed);
      }
        if (GraceOpen(2048) == -1) {
        fprintf(stderr, "Can't run Grace. \n");
        exit(1);
    }
    GracePrintf("s0 on");
    GracePrintf("s0 symbol 1");
    GracePrintf("s0 line linestyle 1");
    GracePrintf("s0 symbol size 0.75");
    GracePrintf("s1 on");
    GracePrintf("s1 symbol 2");
    GracePrintf("s1 line linestyle 1");
    GracePrintf("s1 symbol size 0.75");
    for(unsigned int i=0;i<xi.size();i++){
                GracePrintf("s0 point %g,%g",xi.at(i),yi.at(i));
                GracePrintf("s1 point %g,%g",xi.at(i),xi.at(i)*yi.at(i));
    }
                GracePrintf("autoscale");
                GracePrintf("legend on");
                GracePrintf("legend 0.65,0.45");
                GracePrintf("s0 legend \"Speed per process\"");
                GracePrintf("s1 legend \"Total speed\"");
                GraceFlush();


    return 0;
}





