;;;; cformat.lisp -- convert maxima expressions to c format
;;;; this work in very incomplete as of the writing of this comment

;;;; Copyright (C) 2008 James F. Amundson

;;;; cformat.lisp is free software; you can redistribute it
;;;; and/or modify it under the terms of the GNU General Public
;;;; License as published by the Free Software Foundation; either
;;;; version 2, or (at your option) any later version.

;;;; cformat.lisp is distributed in the hope that it will be
;;;; useful, but WITHOUT ANY WARRANTY; without even the implied
;;;; warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
;;;; See the GNU General Public License for more details.

;;;; You should have received a copy of the GNU General Public License
;;;; along with cformat.lisp; see the file COPYING.  If not,
;;;; write to the Free Software Foundation, Inc., 59 Temple Place -
;;;; Suite 330, Boston, MA 02111-1307, USA.

;;;; Loosely based on fortra.lisp. Copyright statements for fortra.lisp follow:
;;;;  Copyright (c) 1984,1987 by William Schelter,University of Texas
;;;;     All rights reserved
;;;;  (c) Copyright 1980 Massachusetts Institute of Technology

(in-package :maxima)

(macsyma-module cformat)

(defmacro fixp (x) `(typep ,x 'fixnum))
        
(defun cscan (e)
 (cond ((atom e)
                (cond ((fixp e) (float e))
                    ((eq e '$%i) (list '(mprogn) 0.0 1.0)) ;; %I is (0,1)
                    ; ((eq e '$%pi) '$pi)
                    (t e)))
       ((member 'array (cdar e) :test #'eq) e)
       (t (let ((op (caar e)))
            (cond ((eq op '%log) (list '(%alog simp) (cscan (cadr e))))
                  ((eq op 'mexpt)
                   (let ((expon (caddr e)) (mybase (cadr e)))
                     (cond ((eq mybase '$%e) (list '($exp simp) (cscan expon)))
                           ((alike1 expon 1//2) (list '(%sqrt simp) (cscan mybase)))
                           ((alike1 expon -1//2)
                            (list '(mquotient simp) 1 (list '(%sqrt simp) (cscan mybase))))
                           ((and (eq expon 2) (atom mybase))
                                (list '(mtimes simp) mybase mybase))
                           ((and (eq expon 3) (atom mybase))
                                (list '(mtimes simp) mybase mybase mybase))
                           ((and (eq expon -1) (atom mybase))
                                (list '(mquotient simp) 1 mybase))
                           ((and (eq expon -2) (atom mybase))
                                (list '(mquotient simp) 1 (list '(mtimes) mybase mybase)))
                           ((and (eq expon 3) (atom mybase))
                                (list '(mquotient simp) 1 (list '(mtimes) mybase mybase mybase)))
                           ; (t (list (car e)
                                    ; (cscan mybase)
                                    ; (cond ((fixp expon) expon)
                                          ; (t (cscan expon))))))))
                           (t (list '($pow)
                                    (cscan mybase)
                                    (cond ((fixp expon) expon)
                                          (t (cscan expon))))))))
;jfa commented-out becuase we don't have rem-value                  ((eq op 'rat) (rem-value e))
                  ((eq op 'mrat) (cscan (ratdisrep e)))
                  ;;  complex numbers to f77 syntax a+b%i ==> (a,b)
                  ((and (member op '(mtimes mplus) :test #'eq)
                        ((lambda (a) 
                           (and (numberp (cadr a))
                                (numberp (caddr a))
                                (not (zerop1 (cadr a)))
                                (list '(mprogn) (caddr a) (cadr a))))
                         (simplify ($bothcoef e '$%i)))))
                  ((and (eq op 'mtimes) (equal -1 (cadr e)))
                   `((mminus simp) ,(cond ((cdddr e)
                                           (do ((ele (cddr e) (cdr ele))
                                                (nl (list '(mtimes simp))
						    `(,@nl ,(cscan (car ele)))))
                                               ((null ele) nl)))
                                          (t (cscan (caddr e))))))
                  ((and (eq op 'mquotient) (member (cadr e) '(1 -1)))
                   `((mquotient simp) ,(cadr e) ,(cscan (caddr e))))
                  (t (do ((ele (cdr e) (cdr ele))
                          (nl nil `(,@nl ,(cscan (car ele)))))
                         ((null ele) (cons (car e) nl)))))))))

(defun cformat-print (x
		  &aux
		  ;; This is a poor way of saying that array references
		  ;; are to be printed with parens instead of brackets.
		  (lb #\()
		  (rb #\)))
  ;; Restructure the expression for displaying.
  (setq x (cscan x))
  ;; Linearize the expression using MSTRING.  Some global state must be
  ;; modified for MSTRING to generate using Fortran syntax.  This must be
  ;; undone so as not to modifiy the toplevel behavior of MSTRING.
  (unwind-protect
       (defprop mexpt msize-infix grind)
    (defprop mminus 100 lbp)

    (defprop msetq (#\:) strsym)
    (setq x (mstring x))
    ;; Make sure this gets done before exiting this frame.
    (defprop mexpt msz-mexpt grind)
    (remprop 'mminus 'lbp))
  (do ((char 0 (1+ char))
       (line ""))
      ((>= char (length x)))
    (setf line (concatenate 'string line (make-sequence
					  'string 1
					  :initial-element (nth char x))))
    (if (>= (length line) 65)
	(let ((break_point -1))
	  (mapc #'(lambda (x)
		    (let ((p (search x line :from-end t)))
		      (if (and p (> p 0))
			  (setf break_point p))))
		'("+" "-" "*" "/"))
	  (incf break_point)
	  (if (= break_point 0)
	      (progn (princ line) (setf line "     "))
	      (progn
		(princ (subseq line 0 break_point))
		(terpri)
		(setf line (concatenate 'string ""
					(subseq line break_point
						(length line))))))))
    (if (and (= char (1- (length x))) (not (equal line "     ")))
	(princ line)))
  (terpri)
  '$done)

(defun pyformat-print (x
		  &aux
		  ;; This is a poor way of saying that array references
		  ;; are to be printed with parens instead of brackets.
		  (lb #\()
		  (rb #\)))
  ;; Restructure the expression for displaying.
  (setq x (fortscan x))
  ;; Linearize the expression using MSTRING.  Some global state must be
  ;; modified for MSTRING to generate using Fortran syntax.  This must be
  ;; undone so as not to modifiy the toplevel behavior of MSTRING.
  (unwind-protect
       (defprop mexpt msize-infix grind)
    (defprop mminus 100 lbp)

    (defprop msetq (#\:) strsym)
    (setq x (mstring x))
    ;; Make sure this gets done before exiting this frame.
    (defprop mexpt msz-mexpt grind)
    (remprop 'mminus 'lbp))
  (do ((char 0 (1+ char))
       (line ""))
      ((>= char (length x)))
    (setf line (concatenate 'string line (make-sequence
					  'string 1
					  :initial-element (nth char x))))
    (if (>= (length line) 65)
	(let ((break_point -1))
	  (mapc #'(lambda (x)
		    (let ((p (search x line :from-end t)))
		      (if (and p (> p 0))
			  (setf break_point p))))
		'("+" "-" "*" "/"))
	  (incf break_point)
	  (if (= break_point 0)
	      (progn (princ line) (setf line "     "))
	      (progn
		(princ (subseq line 0 break_point))
		(princ " \\")
		(terpri)
		(setf line (concatenate 'string ""
					(subseq line break_point
						(length line))))))))
    (if (and (= char (1- (length x))) (not (equal line "     ")))
	(princ line)))
  (terpri)
  '$done)

(defun octformat-print (x
		  &aux
		  ;; This is a poor way of saying that array references
		  ;; are to be printed with parens instead of brackets.
		  (lb #\()
		  (rb #\)))
  ;; Restructure the expression for displaying.
  (setq x (fortscan x))
  ;; Linearize the expression using MSTRING.  Some global state must be
  ;; modified for MSTRING to generate using Fortran syntax.  This must be
  ;; undone so as not to modifiy the toplevel behavior of MSTRING.
  (unwind-protect
       (defprop mexpt msize-infix grind)
    (defprop mminus 100 lbp)

    (defprop msetq (#\:) strsym)
    (setq x (mstring x))
    ;; Make sure this gets done before exiting this frame.
    (defprop mexpt msz-mexpt grind)
    (remprop 'mminus 'lbp))
  (do ((char 0 (1+ char))
       (line ""))
      ((>= char (length x)))
    (setf line (concatenate 'string line (make-sequence
					  'string 1
					  :initial-element (nth char x))))
    (if (>= (length line) 65)
	(let ((break_point -1))
	  (mapc #'(lambda (x)
		    (let ((p (search x line :from-end t)))
		      (if (and p (> p 0))
			  (setf break_point p))))
		'("+" "-" "*" "/"))
	  (incf break_point)
	  (if (= break_point 0)
	      (progn (princ line) (setf line "     "))
	      (progn
		(princ (subseq line 0 break_point))
		(princ " ...")
		(terpri)
		(setf line (concatenate 'string ""
					(subseq line break_point
						(length line))))))))
    (if (and (= char (1- (length x))) (not (equal line "     ")))
	(princ line)))
  (terpri)
  '$done)
(defmspec $cformat (l)
  (setq l (fexprcheck l))
  (let ((value (strmeval l)))
    (cond ((msetqp l) (setq value `((mequal) ,(cadr l) ,(meval l)))))
    (cond ((and (symbolp l) ($matrixp value))
	   ($fortmx l value))
	  ((and (not (atom value)) (eq (caar value) 'mequal)
		(symbolp (cadr value)) ($matrixp (caddr value)))
	   ($fortmx (cadr value) (caddr value)))
	  (t (cformat-print value)))))

(defmspec $pyformat (l)
  (setq l (fexprcheck l))
  (let ((value (strmeval l)))
    (cond ((msetqp l) (setq value `((mequal) ,(cadr l) ,(meval l)))))
    (cond ((and (symbolp l) ($matrixp value))
	   ($fortmx l value))
	  ((and (not (atom value)) (eq (caar value) 'mequal)
		(symbolp (cadr value)) ($matrixp (caddr value)))
	   ($fortmx (cadr value) (caddr value)))
	  (t (pyformat-print value)))))

(defmspec $octformat (l)
  (setq l (fexprcheck l))
  (let ((value (strmeval l)))
    (cond ((msetqp l) (setq value `((mequal) ,(cadr l) ,(meval l)))))
    (cond ((and (symbolp l) ($matrixp value))
	   ($fortmx l value))
	  ((and (not (atom value)) (eq (caar value) 'mequal)
		(symbolp (cadr value)) ($matrixp (caddr value)))
	   ($fortmx (cadr value) (caddr value)))
	  (t (octformat-print value)))))

(defmspec $jfa (l)
    (format t "jfa: ~a~%" l)
    (setq l (fexprcheck l))
    (format t "jfa2: ~a~%" l)
    (setq l (strmeval l))
    (format t "jfa3: ~a~%" l)
    (setq l (cscan l))
    (format t "jfa4: ~a~%" l)
    (cformat-print l)
    )
    