import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
plt.rcParams["font.family"]="serif"
fig,ax=plt.subplots(figsize=(10,4.2)); ax.set_xlim(0,100); ax.set_ylim(0,42); ax.axis("off")
def box(x,y,w,h,t,fc,fs=8.5,bold=False):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.3,rounding_size=1.2",fc=fc,ec="black",lw=1.1))
    ax.text(x+w/2,y+h/2,t,ha="center",va="center",fontsize=fs,fontweight="bold" if bold else "normal")
def arr(x1,y1,x2,y2,t=None,off=(0,1.4),style="-|>",color="black"):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle=style,mutation_scale=11,lw=1.3 if color!="black" else 1.1,color=color))
    if t: ax.text((x1+x2)/2+off[0],(y1+y2)/2+off[1],t,ha="center",fontsize=7.5,style="italic")
# inputs
box(1,26,14,9,"Workflow DAG\ntasks, files","#eeeeee")
box(1,8,14,9,"VM pool\nprices, MIPS,\nbandwidth","#eeeeee")
# decoder
box(20,6,26,31,"","#dce9f6")
ax.text(33,34,"Simulator-faithful decoder",ha="center",fontsize=9,fontweight="bold")
box(22,24.5,22,6,"FIFO queues, stage-in billed,\nper-VM file cache","white",7.8)
box(22,15.5,22,6.5,"Nominal mode (cheap):\n1 decode","#e4f2e4",7.8)
box(22,8,22,6.5,"Scenario mode (expensive):\nS decodes, mean or CVaR","#f7dede",7.8)
arr(15,30.5,20,30.5); arr(15,12.5,20,12.5)
# search loop
ax.add_patch(FancyBboxPatch((50,3),49,36,boxstyle="round,pad=0.3,rounding_size=1.5",fc="none",ec="#555555",lw=1.4,ls="--"))
ax.text(97.5,4.2,"LIWSA+ search, one generation",ha="right",fontsize=8.5,fontweight="bold")
A,B,C=52,68.5,85
box(C,29,13,8,"1 Parents\n(tournament)","#fbefd3",7.8)
box(B,29,14,8,"2 Draw K x lambda\ncandidates","#fbefd3",7.8)
box(A,29,14,8,"3 Nominal decode,\nrank (ridge optional)","#e4f2e4",7.4)
box(A,14,14,8,"4 Best lambda:\nscenario evaluation","#f7dede",7.6)
box(B,14,14,8,"5 Survival:\nrank + crowding","#e4f2e4",7.8)
box(C,14,13,8,"6 Anytime archive\n(scored at the end)","#e4f2e4",7.2)
arr(C,33,B+14,33); arr(B,33,A+14,33); arr(A+7,29,A+7,22); arr(A+14,18,B,18); arr(B+14,18,C,18); arr(C+6.5,22,C+6.5,29)
ax.text(C+5.8,25.5,"next\ngeneration",fontsize=6.8,style="italic",va="center",ha="right")
arr(44,19.3,A,33,color="#2e7d32")
arr(A,17,44,11.2,color="#c62828")
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import ANALYSIS
(ANALYSIS/"figures").mkdir(exist_ok=True)
fig.savefig(ANALYSIS/"figures"/"fig_framework.pdf",bbox_inches="tight"); fig.savefig(ANALYSIS/"figures"/"fig_framework.png",dpi=150,bbox_inches="tight")
