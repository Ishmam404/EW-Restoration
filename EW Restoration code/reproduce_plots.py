"""Reproduce the two figures in joint_fit_paper_note.pdf from saved matrices."""
from pathlib import Path
import json, os
ROOT=Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'tmp/mpl-cache'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
fits=json.loads((ROOT/'output/joint-sd-fit/fit_results.json').read_text())
TMP=ROOT/'output/pdf'
TMP.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'cmr10','mathtext.fontset':'cm',
    'axes.formatter.use_mathtext':True,'font.size':9,'axes.labelsize':11,
    'axes.linewidth':.7,'pdf.fonttype':42,'savefig.facecolor':'white',
    'axes.unicode_minus':False})
level=2.30; checks=[]

def ellipse(F):
    F=np.asarray(F); t=np.linspace(0,2*np.pi,3001)
    pts=np.sqrt(level)*np.linalg.cholesky(np.linalg.inv(F))@np.array([np.cos(t),np.sin(t)])
    checks.append(float(np.max(np.abs(np.einsum('in,ij,jn->n',pts,F,pts)-level))))
    return .01*pts[1],pts[0]

def axes():
    fig,ax=plt.subplots(figsize=(3.45,3.38))
    fig.subplots_adjust(left=.165,right=.985,bottom=.145,top=.91)
    ax.axhline(0,color='.4',lw=.7,ls=(0,(2,3)))
    ax.axvline(0,color='.8',lw=.5)
    ax.plot(0,0,'o',color='black',ms=3,zorder=8)
    ax.set(xlim=(-.033,.033),ylim=(-10,10),xlabel=r'$C_{Hq}^{(3)}\ [\mathrm{TeV}^{-2}]$',ylabel=r'$n_7/10^{-4}$')
    ax.xaxis.set_major_locator(MultipleLocator(.02))
    ax.xaxis.set_minor_locator(MultipleLocator(.01))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.yaxis.set_major_locator(MultipleLocator(5))
    ax.yaxis.set_minor_locator(MultipleLocator(1))
    ax.tick_params(which='both',direction='in',top=True,right=True)
    return fig,ax

fig,ax=axes(); colors={'S':'#b56435','D':'#209487','joint':'#174b73'}
for name in ['S','D']:
    F=np.asarray(fits[f'{name}_only_2_3000']['F'])
    if name=='S':
        x=np.linspace(-10,10,2001)
        disc=F[1,1]*level-np.linalg.det(F)*x*x
        ys=[(-F[0,1]*x+s*np.sqrt(disc))/F[1,1] for s in [-1,1]]
        ax.fill_betweenx(x,.01*ys[0],.01*ys[1],color=colors[name],alpha=.10)
        for y in ys:
            ax.plot(.01*y,x,color=colors[name],lw=1.2)
            checks.append(float(np.max(np.abs(F[0,0]*x*x+2*F[0,1]*x*y+F[1,1]*y*y-level))))
    else:
        y=np.linspace(-3.3,3.3,2001)
        disc=F[0,0]*level-np.linalg.det(F)*y*y
        xs=[(-F[0,1]*y+s*np.sqrt(disc))/F[0,0] for s in [-1,1]]
        ax.fill_between(.01*y,xs[0],xs[1],color=colors[name],alpha=.10)
        for x in xs:
            ax.plot(.01*y,x,color=colors[name],lw=1.2,ls='--')
            checks.append(float(np.max(np.abs(F[0,0]*x*x+2*F[0,1]*x*y+F[1,1]*y*y-level))))
x,y=ellipse(fits['2_3000_independent']['F'])
ax.fill(x,y,color=colors['joint'],alpha=.13); ax.plot(x,y,color=colors['joint'],lw=1.6)
ax.legend(handles=[Line2D([],[],color=colors['S'],lw=1.2,label=r'$S$ only'),Line2D([],[],color=colors['D'],lw=1.2,ls='--',label=r'$D$ only'),Line2D([],[],color=colors['joint'],lw=1.6,label=r'Joint $S,D$')],loc='upper right',fontsize=8,framealpha=.95,borderpad=.35)
ax.set_title(r'$M_{\max}=3\ \mathrm{TeV}$',fontsize=10,pad=7)

for ext in ['pdf','svg','png']:
    fig.savefig(TMP/f'figure1_sd_comparison.{ext}',dpi=300)
plt.close(fig)
fig,ax=axes()
for cut,color,style in [(500,'#b56435','-'),(1000,'#8062a2','-'),(2000,'#209487','--'),(3000,'#174b73','-')]:
    x,y=ellipse(fits[f'2_{cut}_independent']['F'])
    ax.plot(x,y,color=color,lw=1.4,ls=style,label=rf'${cut/1000:g}\ \mathrm{{TeV}}$')
    if cut==3000: ax.fill(x,y,color=color,alpha=.10)
ax.legend(title=r'$M_{\max}$',loc='upper right',fontsize=8,title_fontsize=8,framealpha=.95,borderpad=.35,labelspacing=.3)
ax.set_title(r'Joint $S,D$ fit',fontsize=10,pad=7)

for ext in ['pdf','svg','png']:
    fig.savefig(TMP/f'figure2_mass_cut_scan.{ext}',dpi=300)
plt.close(fig)
assert max(checks)<1e-10


(TMP/'plot_validation.json').write_text(json.dumps({'max_contour_residual':max(checks),'delta_chi2':level,'source':'output/joint-sd-fit/fit_results.json'},indent=2)+'\n')
print('Both figures saved in',TMP)
