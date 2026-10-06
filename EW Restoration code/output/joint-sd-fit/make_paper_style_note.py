"""Build a one-page plot-focused note from the saved joint-fit results.

All contours are regenerated from fit_results.json. No integration or fit is
changed. Computer Modern typography and vector figures follow the draft style.
"""
from pathlib import Path
import sys, os, json, re, hashlib, logging
logging.getLogger("fontTools").setLevel(logging.ERROR)
from io import BytesIO
ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'tmp/joint-mpl-cache'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import mathtext
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.enums import TA_JUSTIFY, TA_CENTER
from reportlab.platypus import Paragraph, Spacer, Flowable, Frame, Table, TableStyle
from pypdf import PdfReader, PdfWriter, Transformation

TMP=ROOT/'tmp/pdfs/joint-paper-note'; TMP.mkdir(parents=True,exist_ok=True)
OUT=ROOT/'output/pdf/joint_fit_paper_note.pdf'
fits=json.loads((HERE/'fit_results.json').read_text())
fonts=Path(matplotlib.get_data_path())/'fonts/ttf'
for name, fn in [('CM','cmr10.ttf'),('CMB','cmb10.ttf'),('CMI','cmti10.ttf')]:
    pdfmetrics.registerFont(TTFont(name,str(fonts/fn)))
pdfmetrics.registerFontFamily('CM',normal='CM',bold='CMB',italic='CMI',boldItalic='CMB')
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
fig.savefig(TMP/'figure1.pdf');plt.close(fig)
fig,ax=axes()
for cut,color,style in [(500,'#b56435','-'),(1000,'#8062a2','-'),(2000,'#209487','--'),(3000,'#174b73','-')]:
    x,y=ellipse(fits[f'2_{cut}_independent']['F'])
    ax.plot(x,y,color=color,lw=1.4,ls=style,label=rf'${cut/1000:g}\ \mathrm{{TeV}}$')
    if cut==3000: ax.fill(x,y,color=color,alpha=.10)
ax.legend(title=r'$M_{\max}$',loc='upper right',fontsize=8,title_fontsize=8,framealpha=.95,borderpad=.35,labelspacing=.3)
ax.set_title(r'Joint $S,D$ fit',fontsize=10,pad=7)
fig.savefig(TMP/'figure2.pdf');plt.close(fig)
assert max(checks)<1e-10

# Text is Computer Modern; inline equations are rendered at 600 dpi and display
# equations/figures are embedded as vector PDF overlays.
parser=mathtext.MathTextParser('path'); overlays=[]; equation_count=0

def math_asset(expr,size=10,kind='png'):
    key=hashlib.sha1(f'{size}:{kind}:{expr}'.encode()).hexdigest()[:16]
    path=TMP/f'math-{key}.{kind}'
    prop=FontProperties(size=size)
    parsed=parser.parse('$'+expr+'$',dpi=72,prop=prop)
    if not path.exists():
        mathtext.math_to_image('$'+expr+'$',str(path),prop=prop,dpi=600,format=kind,color='black')
    if kind=='pdf':
        box=PdfReader(path).pages[0].mediabox
        return path,float(box.width),float(box.height),parsed.depth
    # math_to_image uses ceil(width), ceil(height).
    return path,parsed.width,parsed.height,parsed.depth

def rich(text,size=10):
    def repl(match):
        path,w,h,depth=math_asset(match.group(1),size)
        return f'<img src="{path}" width="{w}" height="{h}" valign="{-depth}"/>'
    return re.sub(r'\$(.+?)\$',repl,text)

styles={
 'body':ParagraphStyle('body',fontName='CM',fontSize=10,leading=12.2,alignment=TA_JUSTIFY,firstLineIndent=10,spaceAfter=4),
 'caption':ParagraphStyle('caption',fontName='CM',fontSize=9,leading=10.7,alignment=TA_JUSTIFY,spaceAfter=6),
 'abstract':ParagraphStyle('abstract',fontName='CM',fontSize=9.3,leading=11.4,alignment=TA_JUSTIFY),
 'head':ParagraphStyle('head',fontName='CMB',fontSize=9.5,leading=12,alignment=TA_CENTER,spaceBefore=8,spaceAfter=9),
 'small':ParagraphStyle('small',fontName='CM',fontSize=9,leading=10.8,alignment=TA_JUSTIFY,spaceAfter=4),
}

def P(text,kind='body'):
    return Paragraph(rich(text,styles[kind].fontSize),styles[kind])

def H(text): return Paragraph(text,styles['head'])

class Vector(Flowable):
    def __init__(self,path,width,height=None,label=None):
        super().__init__();self.path=path;self.width=width;self.label=label
        page=PdfReader(path).pages[0]; self.original=(float(page.mediabox.width),float(page.mediabox.height))
        self.height=height or width*self.original[1]/self.original[0]
    def draw(self):
        x,y=self.canv.absolutePosition(0,0)
        overlays.append((self.canv.getPageNumber()-1,self.path,x,y,self.width,self.height))
        if self.label:
            self.canv.setFont('CM',10);self.canv.drawRightString(243,self.height/2-3,self.label)

class Equation(Vector):
    def __init__(self,expr):
        global equation_count
        equation_count+=1
        path,w,h,depth=math_asset(expr,10,'pdf')
        super().__init__(path,min(w,216),label=f'({equation_count})')
        self._spaceBefore=4;self._spaceAfter=7
    def draw(self):
        x,y=self.canv.absolutePosition((243-self.width)/2-5,0)
        overlays.append((self.canv.getPageNumber()-1,self.path,x,y,self.width,self.height))
        self.canv.setFont('CM',10);self.canv.drawRightString(243,self.height/2-3,self.label)
    def getSpaceBefore(self):return self._spaceBefore
    def getSpaceAfter(self):return self._spaceAfter

base=TMP/'base.pdf'; c=canvas.Canvas(str(base),pagesize=(612,792))
c.setTitle('Joint constraints on common and decorrelating electroweak directions')
c.setAuthor('');c.setSubject('S-only, D-only and joint fits; invariant-mass-cut dependence')

def frame(items,x,y,w,h,name):
    remaining=list(items)
    f=Frame(x,y,w,h,leftPadding=0,rightPadding=0,topPadding=0,bottomPadding=0,id=name)
    f.addFromList(remaining,c)
    if remaining:
        raise RuntimeError(f'{name}: {len(remaining)} flowables did not fit; first={getattr(remaining[0],"text",type(remaining[0]).__name__)}')

# One-page note: retain the plots and focus the prose on their interpretation.
c.setFont('CMB',14)
c.drawCentredString(306,750,'Joint fit of common and decorrelating directions')
intro=(r'The plots show expected constraints on $C_{Hq}^{(3)}$ and $n_7$ from normalized longitudinal '
 r'$WZ$ and $Wh$ rates. The observables $S=(\mu_Z+\mu_h)/2-1$ and $D=(\mu_Z-\mu_h)/2$ '
 r'measure their average shift and difference, respectively; $\mu$ denotes the rate divided by its SM prediction.')
frame([P(intro)],54,661,504,65,'intro')
frame([Vector(TMP/'figure1.pdf',243)],54,412,243,240,'figure1')
frame([Vector(TMP/'figure2.pdf',243)],315,412,243,240,'figure2')
cap1=(r'FIG. 1. $S$-only (orange), $D$-only (green dashed), and joint $S,D$ (blue) constraints '
 r'for $M<M_{\max}=3\ \mathrm{TeV}$. Both coefficients float. The single-observable bands extend beyond the frame.')
cap2=(r'FIG. 2. Joint $S,D$ constraints for $M_{\max}=0.5,1,2,3\ \mathrm{TeV}$. '
 r'Each contour uses a separate cumulative mass-selected sample, with fixed transverse-momentum bins and fractional errors.')
frame([P(cap1,'caption')],54,337,243,68,'caption1')
frame([P(cap2,'caption')],315,337,243,68,'caption2')
body=[
P(r'In Fig. 1, the nearly vertical $S$-only band mainly constrains $C_{Hq}^{(3)}$, while the nearly horizontal $D$-only band mainly constrains $n_7$. Combining them, including their covariance, produces a compact, tilted ellipse: the two observables constrain complementary directions.'),
P(r'In Fig. 2, the contours shrink as the upper mass cut increases. The largest improvement occurs between 0.5 and 2 TeV; the 2 and 3 TeV contours nearly overlap. From 0.5 to 3 TeV, the profiled uncertainties improve by factors of about 1.9 for $n_7$ and 1.8 for $C_{Hq}^{(3)}$, with only a further 2% improvement from 2 to 3 TeV.'),
P(r'All boundaries correspond to $\Delta\chi^2=2.30$ (68.3% joint confidence for two parameters). The black point marks the SM; the horizontal dotted line, $n_7=0$, is the selected SMEFT direction. Every contour is centered on the SM, so these are expected sensitivities rather than evidence for a nonzero deformation.'),
P(r'Benchmark: $\sqrt{s}=14\ \mathrm{TeV}$, nCTEQ15, interference-order templates, and $p_T$ bins of 100-200 and 200-1500 GeV. Fractional channel errors are $\epsilon_Z=(0.088,0.091)$ and $\epsilon_h=(0.15,0.15)$, with uncorrelated channel and bin measurements. The mass scan holds these errors fixed.','small')]
frame(body,54,55,504,270,'description')
c.save()
reader=PdfReader(base)
for page_index,path,x,y,w,h in overlays:
    src=PdfReader(path).pages[0]
    sx=w/float(src.mediabox.width);sy=h/float(src.mediabox.height)
    reader.pages[page_index].merge_transformed_page(src,Transformation().scale(sx,sy).translate(x,y))
writer=PdfWriter()
for page in reader.pages:writer.add_page(page)
writer.add_metadata({'/Title':'Joint constraints on common and decorrelating electroweak directions','/Subject':'Two-figure note: S-only, D-only, joint fits and upper invariant-mass scan'})
with open(OUT,'wb') as f:writer.write(f)
assert len(PdfReader(OUT).pages)==1
(HERE/'paper_note_validation.json').write_text(json.dumps({'pages':1,'max_contour_residual':max(checks),'source':'fit_results.json','overlays':len(overlays),'delta_chi2':level},indent=2)+'\n')
print('Created',OUT,'with',len(overlays),'vector overlays; contour residual',max(checks))
