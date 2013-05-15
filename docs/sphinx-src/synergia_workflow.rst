synergia_workflow module
========================

The synergia_workflow module provides support for parameterizing simulations, 
interacting with batch systems and keeping track of simulation runs.

-------------------
Python-only Classes
-------------------

.. automodule:: synergia_workflow
   
.. autoclass:: synergia_workflow.Options
   :members:

.. autoclass:: synergia_workflow.Job_manager
   :members:

.. autoclass:: synergia_workflow.Override

Built-in Job_manager options
----------------------------

=========================== ============= =================================== ====================
name                        type          default                             description         
=========================== ============= =================================== ====================
account                     :code:`str`   None                                Batch system account, default=None
createjob                   :code:`bool`  False                               Whether to create new job directory, default=False
jobdir                      :code:`str`   run                                 Job directory, default=run
multitemplate               :code:`str`   multijob                            Base filename for multijob templates to look for in search path, default=multijob
multitemplatepath           :code:`str`   None                                Full path to the base name for multijob templates, default=None
numproc                     :code:`int`   1                                   Number of processes, default=1
overwrite                   :code:`bool`  False                               Whether to overwrite existing job directory, default=False
procspernode                :code:`int`   1                                   Number of processes per node, default=1
queue                       :code:`str`   None                                Batch system queue, default=None
resumemultitemplate         :code:`str`   resumemultijob                      Base filename for resume multijob templates to look for in search path, default=resumemultijob
resumemultitemplatepath     :code:`str`   None                                Full path to the base name for resume multijob templates, default=None
resumetemplate              :code:`str`   resumejob                           Filename for resume job template to look for in search path, default=resumejob
resumetemplatepath          :code:`str`   None                                Full path to a resume job template, default=None
run                         :code:`bool`  False                               Whether to immediately run job, default=False
setupsh                     :code:`str`   /usr/local/share/synergia/setup.sh  Path to Synergia2 setup.sh file, default=/usr/local/share/synergia/setup.sh
submit                      :code:`bool`  False                               Whether to immediately submit job, default=False
synergia_executable         :code:`str`   synergia                            Name or path of synergia executable, default=synergia
synergia_resume_executable  :code:`str`   synergia-pyresume                   Name or path of synergia resume executable, default=synergia-pyresume
template                    :code:`str`   job                                 Filename for job template to look for in search path, default=job
templatepath                :code:`str`   None                                Full path to a job template, default=None
walltime                    :code:`str`   None                                Limit job to given wall time, default=None
=========================== ============= =================================== ====================

