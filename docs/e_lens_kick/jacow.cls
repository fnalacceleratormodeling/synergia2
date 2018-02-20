%%
%% This file has been developed as a common template for papers
%% destined for electronic production using Adobe Acrobat
%% software for Accelerator Conferences
%%
%% See the JACo Website for more information
%%
%%       http://jacow.org/
%%
%  v0.1-1.3
%  Special thanks to John Jowett and Michel Goossens from CERN and
%  Martin Comyn at TRIUMF for their significant contributions to
%  this class file over the period 1996 to 2000.
%                                                 John Poole
%                                                 March 2000
%  v1.4
%  JAC2001.cls is a modified version of JAC2000.cls to produce indented
%  first paragraphs after section, subsection and subsubsection headings.
%                                                 Martin Comyn  April 2001
%
%  v1.5
%  This is a modified version of JAC2003.cls to adjust space around
%  section and subsection headers to be more consistent with JACOW Word
%  templates.                                     Todd Satogata  March 2011
%
% v 1.6
% This is a complectly rewritten version of JAC2003.cls which needs a current
% TeX-System to run.
%                                                 Ulrike Fischer, November 2013
%
% v 1.7
% - small change to correct the text block inside JACoW's magic red borders for
%   a4paper (aca4); top has been set 18.5mm (19mm is defined in the template but
%   leaves descenders outside the lower y margin).
% - duplicate {boxit} removed
%                                                 Volker RW Schaa, 14 April 2014
%
% v1.8
% - added setup for \micro sign which disappears when using XeTeX or LuaTeX
%   with unicode-math.                             Ulrike Fischer, 21 April 2014
%
% v1.9
% - fixed the pdfLaTeX warnings for the text/math-micro hack
%                                                  Ulrike Fischer, 22 April 2014
%
% v1.91
% - Ligatures=TeX switch introduced to accommodate
%                                                  Ulrike Fischer, 22 April 2014
%
% v1.92
% - settings for top margin have to be different in A4 and letter to accommodate
%   JACoW's PitStop Action List. This was found after receiving Plamen Hopchev's
%   email about margins and testing the workflow with cropping the bounding box
%   which starts at the lower left edge and not at the top (see graphic JACpic_mc
%   in the template for measures).
%                                                 Volker RW Schaa, 29 April 2014
% v1.93
% - setting the bottom margin (19mm) without top solves the problem for different
%   A4/Letter settings. This was already the default in v1.6. Pointed out by
%   Plamen Hopchev. To accommodate the descenders the bottom margin has been set
%   to 56pt now.
%                                                   Volker RW Schaa, 01 May 2014
%
% v1.94
% - the micro sign in UTF-8 prevents ASCII format of the cls file. Ulrike pointed
%   out a hack in http://tex.stackexchange.com/questions/172968/hide-notation-from-pdftex
%   which is now introduced.
%                                                   Volker RW Schaa, 02 May 2014
%
\def\fileversion{1.94}
\def\filedate{2014/05/02}
\def\docdate {2014/05/02}

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{jacow}[\filedate\space Version \fileversion]

\typeout{------------------------------------------------------------------------}
\typeout{LaTeX2e Class file for Accelerator Conference publication for LaTeX2e users}
\typeout{ }
\typeout{Use the boxit option to draw a box on page showing the correct margins}
\typeout{ }
\typeout{Itemize, Enumerate and Description environments are compact versions}
\typeout{------------------------------------------------------------------------}
\typeout{ }

\newif\ifjacowbiblatex
\newif\ifjacowrefpage

\DeclareOption{acus}{%
   \PassOptionsToPackage{paper=letterpaper}{geometry}
   \typeout{Setup for US LETTER PAPER}}

\DeclareOption{letterpaper}{%
   \PassOptionsToPackage{paper=letterpaper}{geometry}
   \typeout{Setup for US LETTER PAPER}}

\DeclareOption{a4paper}{%
    \PassOptionsToPackage{paper=a4paper}{geometry}
    \typeout{Setup for A4 PAPER}}

\DeclareOption{aca4}{%
    \PassOptionsToPackage{paper=a4paper}{geometry}
    \typeout{Setup for A4 PAPER}}

\DeclareOption{boxit}{\PassOptionsToPackage{showframe}{geometry}}

\DeclareOption{biblatex}{\jacowbiblatextrue}

\DeclareOption{refpage}{\jacowrefpagetrue}

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}

\ExecuteOptions{aca4}
\ProcessOptions
\RequirePackage{fix-cm}
\LoadClass[10pt,twocolumn]{article}
\RequirePackage{flushend}
% Tools:
\RequirePackage{etoolbox}
\RequirePackage{ifxetex}
\RequirePackage{ifluatex}
\RequirePackage{textcase}%

%Add thanks to the list of "\@nonchangecase"-commands from textcase:
\def\@uclcnotmath#1#2#3#4{\begingroup
      #1%
      \def\({$}\let\)\(%
      \def\NoCaseChange##1{\noexpand\NoCaseChange{\noexpand##1}}%
      \@nonchangecase\label
      \@nonchangecase\ref
      \@nonchangecase\ensuremath
      \@nonchangecase\thanks %new
      \@nonchangecase\si %new
      \@nonchangecase\SI %new
      \def\cite##1##{\toks@{\noexpand\cite##1}\@citex}%
      \def\@citex##1{\NoCaseChange{\the\toks@{##1}}}%
      \def\reserved@a##1##2{\let#2\reserved@a}%
      \expandafter\reserved@a\@uclclist\reserved@b{\reserved@b\@gobble}%
      \protected@edef\reserved@a{\endgroup
          \noexpand\@skipmath#3#4$\valign$}%
      \reserved@a}

\RequirePackage[detect-mode,detect-weight]{siunitx}
\RequirePackage{graphicx}
\RequirePackage{booktabs}
\RequirePackage[figureposition=bottom,tableposition=top,skip=5pt]{caption}
\RequirePackage{xcolor}
\RequirePackage{url}
\AtEndPreamble{\RequirePackage[autostyle]{csquotes}}


% Page layout:
\RequirePackage[%
  textwidth=170mm,
  textheight=241mm,
  heightrounded,
  left=20mm,
  bottom=56pt,
  columnsep=5mm,
  noheadfoot,
  nomarginpar,
  twocolumn]
  {geometry}

\columnseprule 0pt
\usepackage[hang]{footmisc}
\setlength{\footnotemargin}{.6em}


\pagestyle{empty}

\ifboolexpr{bool{xetex} or bool{luatex}}
 {}
 { \catcode`\^^^=9
 }
 
\ifboolexpr{bool{xetex} or bool{luatex}}
 { \let\ori@vdots\vdots
   \RequirePackage{unicode-math}
   \AtBeginDocument{\let\vdots\ori@vdots}
   \setmainfont[Ligatures=TeX]{TeX Gyre Termes}
   \setmathfont{TeX Gyre Termes Math}
   \sisetup{
     math-micro = \text{^^^^03bc},
     text-micro = ^^^^03bc
      }
 }
 {
  % Fonts: Times clones
  \RequirePackage{textcomp}
  \RequirePackage[T1]{fontenc}
  \RequirePackage{lmodern}
  \RequirePackage{tgtermes}
  \RequirePackage{newtxmath}
 }

\RequirePackage{microtype}
%Lists

\RequirePackage{enumitem}
\newenvironment{Enumerate}{\begin{enumerate}[nosep]}{\end{enumerate}}
\newenvironment{Itemize}{\begin{itemize}[nosep]}{\end{itemize}}
\newenvironment{Description}{\begin{description}[nosep]}{\end{description}}


%Floatparameter:
\renewcommand{\topfraction}{.95}
\renewcommand{\bottomfraction}{.95}
\renewcommand{\textfraction}{0.1}
\renewcommand{\floatpagefraction}{0.8}


%headings:
% section: Uppercase only for text
\renewcommand{\section}
   {%
    \@startsection{section}{1}{0mm}
       {2.0ex plus 0.8ex minus .1ex}{1.0ex plus .2ex}
       {\normalfont\large\bfseries\mathversion{bold}\centering\MakeTextUppercase}%
   }%

\renewcommand\subsection
  {%
   \@startsection{subsection}{2}{\z@}
    {1.4ex plus .8ex minus .17ex}{0.8ex plus .17ex}
    {\normalfont\large\itshape}%
   }

\renewcommand\subsubsection
 {%
  \@startsection{subsubsection}{3}{\parindent}
  {2.5ex plus .7ex minus .17ex}{-1em}
  {\normalfont\normalsize\bfseries}%
 }

\renewcommand\paragraph
 {%
  \@startsection{paragraph}{4}{\z@}
  {2.5ex plus .7ex minus .17ex}{-1em}
  {\normalfont\normalsize\itshape}%
 }

\renewcommand\subparagraph
 {%
  \@startsection{subparagraph}{4}{\parindent}
  {2.25ex plus .7ex minus .17ex}{-1em}
  {\normalfont\normalsize\bfseries}%
 }

\setcounter{secnumdepth}{0}

% This definition of \maketitle taken from article.sty, and has been
% somewhat modified.

\def\maketitle{\par
 \begingroup
   \def\thefootnote{\fnsymbol{footnote}}
   \def\@makefnmark{\hbox
       to 5pt{$^{\@thefnmark}$\hss}}
   \twocolumn[\@maketitle]
   \@thanks
 \endgroup
 \enlargethispage{\jac@copyrightspace}%
 \setcounter{footnote}{0}
 \let\maketitle\relax
 \let\@maketitle\relax
 \gdef\@thanks{}\gdef\@author{}\gdef\@title{}\let\thanks\relax}

\newlength{\titleblockheight}       % so user can change it if need be
\setlength{\titleblockheight}{3.5cm}

\newlength\titleblockstartskip
\setlength\titleblockstartskip{3pt}


\newlength\titleblockmiddleskip
\setlength\titleblockmiddleskip{1em}

\newlength\titleblockendskip
\setlength\titleblockendskip{1em}


\def\@maketitle{%
  \vskip \titleblockstartskip \centering
  {\Large\bfseries \MakeTextUppercase{\@title} \par}
  \vskip \titleblockmiddleskip               % Vertical space after title.
  {\large\begin{tabular}[t]{@{}c@{}}\@author \end{tabular}\par}
  \vskip \titleblockendskip}


% The \copyrightspace command is used to produce a blank space in the first
% column where a copyright notice may go.  It works by producing
% with \enlargethispage and is inserted by \maketitle.
% The command should be issued in the preamble.

\newcommand\jac@copyrightspace{0pt}
\newcommand\copyrightspace[1][1cm]{\renewcommand\jac@copyrightspace{-#1}}

\ifboolexpr{bool{@titlepage}}
{\renewenvironment{abstract}
 {\list{}{%
    \setlength{\leftmargin}{\dimexpr\textwidth/2-0.75\columnwidth}%
    \setlength{\rightmargin}{\dimexpr-0.75\columnwidth-\columnsep}%
    \setlength{\listparindent}{\parindent}%
    \setlength{\itemsep}{\parskip}%
    \setlength{\itemindent}{\z@}%
    \setlength{\topsep}{\z@}%
    \setlength{\parsep}{\parskip}%
    \setlength{\partopsep}{\z@}%
    \let\makelabel\@gobble
    \setlength{\labelwidth}{\z@}%
    \advance\@listdepth\m@ne   }%
   \item\relax\subsection*{Abstract}}
 {\endlist\clearpage}
}
{%
 \renewenvironment{abstract}
  {\subsection*{Abstract}}
  {\par}
}
\ifboolexpr{bool{jacowbiblatex}}
 {\RequirePackage[style=numeric-comp,sorting=none]{biblatex}
  \renewbibmacro*{url+urldate}{%
    \iffieldundef{url}
     {}
     {\printfield{url}%
      \nopunct}}%
  \DeclareFieldFormat{url}{\url{#1}}
  \DeclareFieldFormat{eprint}{#1}
  \setlength\bibitemsep{0pt}
  \setlength\bibparsep{0pt}
  \setlength\biblabelsep{5pt}
  \ifjacowrefpage\preto\blx@bibliography{\clearpage}\fi
  \AtBeginBibliography{\small\clubpenalty4000\widowpenalty4000}%
 }
 {\RequirePackage{cite}
  % Redefine to use smaller fonts
  \def\thebibliography#1{\setlength{\itemsep}{0pt}\setlength{\parsep}{0pt}%
  \ifjacowrefpage\clearpage\fi
  \section*{REFERENCES\@mkboth
  {REFERENCES}{REFERENCES}}\small\list
  {[\arabic{enumi}]}{\settowidth\labelwidth{[#1]}\leftmargin\labelwidth
    \advance\leftmargin\labelsep
    \usecounter{enumi}}
    \def\newblock{\hskip .11em plus .33em minus .07em}
    \sloppy\clubpenalty4000\widowpenalty4000
    \sfcode`\.=1000\relax}
   \let\endthebibliography=\endlist
  }


%\sloppy
\clubpenalty10000\widowpenalty10000
\flushbottom
%-----------------------------------------------------------------------

%avoid bug of fixltx2e:
%http://www.latex-project.org/cgi-bin/ltxbugs2html?pr=latex/4023
\RequirePackage{fixltx2e}%
\def\@outputdblcol{%
  \if@firstcolumn
    \global\@firstcolumnfalse
    \global\setbox\@leftcolumn\copy\@outputbox
    \splitmaxdepth\maxdimen
    \vbadness\maxdimen
    \setbox\@outputbox\vbox{\unvbox\@outputbox\unskip}%new
    \setbox\@outputbox\vsplit\@outputbox to\maxdimen
    \toks@\expandafter{\topmark}%
    \xdef\@firstcoltopmark{\the\toks@}%
    \toks@\expandafter{\splitfirstmark}%
    \xdef\@firstcolfirstmark{\the\toks@}%
    \ifx\@firstcolfirstmark\@empty
      \global\let\@setmarks\relax
    \else
      \gdef\@setmarks{%
        \let\firstmark\@firstcolfirstmark
        \let\topmark\@firstcoltopmark}%
    \fi
  \else
    \global\@firstcolumntrue
    \setbox\@outputbox\vbox{%
     \hb@xt@\textwidth{%
        \hb@xt@\columnwidth{\box\@leftcolumn \hss}%
        \hfil
        \vrule \@width\columnseprule
        \hfil
       \hb@xt@\columnwidth{\box\@outputbox \hss}}}%
  \@combinedblfloats
    \@setmarks
    \@outputpage
    \begingroup
      \@dblfloatplacement
      \@startdblcolumn
      \@whilesw\if@fcolmade \fi{\@outputpage\@startdblcolumn}%
    \endgroup
  \fi}

\endinput
