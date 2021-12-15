var data = {lines:[
{"lineNum":"    1","line":"/*"},
{"lineNum":"    2","line":"//@HEADER"},
{"lineNum":"    3","line":"// ************************************************************************"},
{"lineNum":"    4","line":"//"},
{"lineNum":"    5","line":"//                        Kokkos v. 3.0"},
{"lineNum":"    6","line":"//       Copyright (2020) National Technology & Engineering"},
{"lineNum":"    7","line":"//               Solutions of Sandia, LLC (NTESS)."},
{"lineNum":"    8","line":"//"},
{"lineNum":"    9","line":"// Under the terms of Contract DE-NA0003525 with NTESS,"},
{"lineNum":"   10","line":"// the U.S. Government retains certain rights in this software."},
{"lineNum":"   11","line":"//"},
{"lineNum":"   12","line":"// Redistribution and use in source and binary forms, with or without"},
{"lineNum":"   13","line":"// modification, are permitted provided that the following conditions are"},
{"lineNum":"   14","line":"// met:"},
{"lineNum":"   15","line":"//"},
{"lineNum":"   16","line":"// 1. Redistributions of source code must retain the above copyright"},
{"lineNum":"   17","line":"// notice, this list of conditions and the following disclaimer."},
{"lineNum":"   18","line":"//"},
{"lineNum":"   19","line":"// 2. Redistributions in binary form must reproduce the above copyright"},
{"lineNum":"   20","line":"// notice, this list of conditions and the following disclaimer in the"},
{"lineNum":"   21","line":"// documentation and/or other materials provided with the distribution."},
{"lineNum":"   22","line":"//"},
{"lineNum":"   23","line":"// 3. Neither the name of the Corporation nor the names of the"},
{"lineNum":"   24","line":"// contributors may be used to endorse or promote products derived from"},
{"lineNum":"   25","line":"// this software without specific prior written permission."},
{"lineNum":"   26","line":"//"},
{"lineNum":"   27","line":"// THIS SOFTWARE IS PROVIDED BY NTESS \"AS IS\" AND ANY"},
{"lineNum":"   28","line":"// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE"},
{"lineNum":"   29","line":"// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR"},
{"lineNum":"   30","line":"// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE"},
{"lineNum":"   31","line":"// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,"},
{"lineNum":"   32","line":"// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,"},
{"lineNum":"   33","line":"// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR"},
{"lineNum":"   34","line":"// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF"},
{"lineNum":"   35","line":"// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING"},
{"lineNum":"   36","line":"// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS"},
{"lineNum":"   37","line":"// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."},
{"lineNum":"   38","line":"//"},
{"lineNum":"   39","line":"// Questions? Contact Christian R. Trott (crtrott@sandia.gov)"},
{"lineNum":"   40","line":"//"},
{"lineNum":"   41","line":"// ************************************************************************"},
{"lineNum":"   42","line":"//@HEADER"},
{"lineNum":"   43","line":"*/"},
{"lineNum":"   44","line":""},
{"lineNum":"   45","line":"// Experimental unified task-data parallel manycore LDRD"},
{"lineNum":"   46","line":""},
{"lineNum":"   47","line":"#ifndef KOKKOS_IMPL_TASKQUEUE_HPP"},
{"lineNum":"   48","line":"#define KOKKOS_IMPL_TASKQUEUE_HPP"},
{"lineNum":"   49","line":""},
{"lineNum":"   50","line":"#include <Kokkos_Macros.hpp>"},
{"lineNum":"   51","line":"#if defined(KOKKOS_ENABLE_TASKDAG)"},
{"lineNum":"   52","line":""},
{"lineNum":"   53","line":"#include <Kokkos_TaskScheduler_fwd.hpp>"},
{"lineNum":"   54","line":"#include <Kokkos_Core_fwd.hpp>"},
{"lineNum":"   55","line":""},
{"lineNum":"   56","line":"#include <Kokkos_MemoryPool.hpp>"},
{"lineNum":"   57","line":""},
{"lineNum":"   58","line":"#include <impl/Kokkos_TaskBase.hpp>"},
{"lineNum":"   59","line":"#include <impl/Kokkos_TaskResult.hpp>"},
{"lineNum":"   60","line":""},
{"lineNum":"   61","line":"#include <impl/Kokkos_Memory_Fence.hpp>"},
{"lineNum":"   62","line":"#include <impl/Kokkos_Atomic_Increment.hpp>"},
{"lineNum":"   63","line":"#include <impl/Kokkos_OptionalRef.hpp>"},
{"lineNum":"   64","line":"#include <impl/Kokkos_LIFO.hpp>"},
{"lineNum":"   65","line":""},
{"lineNum":"   66","line":"#include <string>"},
{"lineNum":"   67","line":"#include <typeinfo>"},
{"lineNum":"   68","line":"#include <stdexcept>"},
{"lineNum":"   69","line":""},
{"lineNum":"   70","line":"//----------------------------------------------------------------------------"},
{"lineNum":"   71","line":"//----------------------------------------------------------------------------"},
{"lineNum":"   72","line":""},
{"lineNum":"   73","line":"namespace Kokkos {"},
{"lineNum":"   74","line":"namespace Impl {"},
{"lineNum":"   75","line":""},
{"lineNum":"   76","line":"/** \\brief  Manage task allocation, deallocation, and scheduling."},
{"lineNum":"   77","line":" *"},
{"lineNum":"   78","line":" *  Task execution is deferred to the TaskQueueSpecialization."},
{"lineNum":"   79","line":" *  All other aspects of task management have shared implementation."},
{"lineNum":"   80","line":" */"},
{"lineNum":"   81","line":"template <typename ExecSpace, typename MemorySpace>"},
{"lineNum":"   82","line":"class TaskQueue : public TaskQueueBase {"},
{"lineNum":"   83","line":" protected:"},
{"lineNum":"   84","line":"  template <class>"},
{"lineNum":"   85","line":"  friend struct TaskQueueSpecialization;"},
{"lineNum":"   86","line":"  template <class, class>"},
{"lineNum":"   87","line":"  friend class TaskQueueSpecializationConstrained;"},
{"lineNum":"   88","line":"  template <class, class>"},
{"lineNum":"   89","line":"  friend class Kokkos::BasicTaskScheduler;"},
{"lineNum":"   90","line":""},
{"lineNum":"   91","line":"  using execution_space = ExecSpace;"},
{"lineNum":"   92","line":"  using memory_space    = MemorySpace;"},
{"lineNum":"   93","line":"  using device_type     = Kokkos::Device<execution_space, memory_space>;"},
{"lineNum":"   94","line":"  using memory_pool     = Kokkos::MemoryPool<device_type>;"},
{"lineNum":"   95","line":"  using task_root_type  = Kokkos::Impl::TaskBase;"},
{"lineNum":"   96","line":"  using team_queue_type = TaskQueue;"},
{"lineNum":"   97","line":""},
{"lineNum":"   98","line":"  struct Destroy {"},
{"lineNum":"   99","line":"    TaskQueue* m_queue;"},
{"lineNum":"  100","line":"    void destroy_shared_allocation();"},
{"lineNum":"  101","line":"  };"},
{"lineNum":"  102","line":""},
{"lineNum":"  103","line":"  //----------------------------------------"},
{"lineNum":"  104","line":""},
{"lineNum":"  105","line":"  enum : int { NumQueue = 3 };"},
{"lineNum":"  106","line":""},
{"lineNum":"  107","line":"  // Queue is organized as [ priority ][ type ]"},
{"lineNum":"  108","line":""},
{"lineNum":"  109","line":"  memory_pool m_memory;"},
{"lineNum":"  110","line":"  task_root_type* volatile m_ready[NumQueue][2];"},
{"lineNum":"  111","line":"  // long                      m_accum_alloc ; // Accumulated number of"},
{"lineNum":"  112","line":"  // allocations"},
{"lineNum":"  113","line":"  int m_count_alloc = 0;  // Current number of allocations","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  114","line":"  int m_max_alloc;        // Maximum number of allocations"},
{"lineNum":"  115","line":"  int m_ready_count;      // Number of ready or executing"},
{"lineNum":"  116","line":""},
{"lineNum":"  117","line":"  //----------------------------------------"},
{"lineNum":"  118","line":""},
{"lineNum":"  119","line":"  ~TaskQueue();"},
{"lineNum":"  120","line":"  TaskQueue()                 = delete;"},
{"lineNum":"  121","line":"  TaskQueue(TaskQueue&&)      = delete;"},
{"lineNum":"  122","line":"  TaskQueue(TaskQueue const&) = delete;"},
{"lineNum":"  123","line":"  TaskQueue& operator=(TaskQueue&&) = delete;"},
{"lineNum":"  124","line":"  TaskQueue& operator=(TaskQueue const&) = delete;"},
{"lineNum":"  125","line":""},
{"lineNum":"  126","line":"  TaskQueue(const memory_pool& arg_memory_pool);"},
{"lineNum":"  127","line":""},
{"lineNum":"  128","line":"  // Schedule a task"},
{"lineNum":"  129","line":"  //   Precondition:"},
{"lineNum":"  130","line":"  //     task is not executing"},
{"lineNum":"  131","line":"  //     task->m_next is the dependence or zero"},
{"lineNum":"  132","line":"  //   Postcondition:"},
{"lineNum":"  133","line":"  //     task->m_next is linked list membership"},
{"lineNum":"  134","line":"  KOKKOS_FUNCTION void schedule_runnable(task_root_type*);"},
{"lineNum":"  135","line":"  KOKKOS_FUNCTION void schedule_aggregate(task_root_type*);"},
{"lineNum":"  136","line":""},
{"lineNum":"  137","line":"  // Reschedule a task"},
{"lineNum":"  138","line":"  //   Precondition:"},
{"lineNum":"  139","line":"  //     task is in Executing state"},
{"lineNum":"  140","line":"  //     task->m_next == LockTag"},
{"lineNum":"  141","line":"  //   Postcondition:"},
{"lineNum":"  142","line":"  //     task is in Executing-Respawn state"},
{"lineNum":"  143","line":"  //     task->m_next == 0 (no dependence)"},
{"lineNum":"  144","line":"  KOKKOS_FUNCTION"},
{"lineNum":"  145","line":"  void reschedule(task_root_type*);"},
{"lineNum":"  146","line":""},
{"lineNum":"  147","line":"  // Complete a task"},
{"lineNum":"  148","line":"  //   Precondition:"},
{"lineNum":"  149","line":"  //     task is not executing"},
{"lineNum":"  150","line":"  //     task->m_next == LockTag  =>  task is complete"},
{"lineNum":"  151","line":"  //     task->m_next != LockTag  =>  task is respawn"},
{"lineNum":"  152","line":"  //   Postcondition:"},
{"lineNum":"  153","line":"  //     task->m_wait == LockTag  =>  task is complete"},
{"lineNum":"  154","line":"  //     task->m_wait != LockTag  =>  task is waiting"},
{"lineNum":"  155","line":"  KOKKOS_FUNCTION"},
{"lineNum":"  156","line":"  void complete(task_root_type*);"},
{"lineNum":"  157","line":""},
{"lineNum":"  158","line":"  KOKKOS_FUNCTION"},
{"lineNum":"  159","line":"  static bool push_task(task_root_type* volatile* const, task_root_type* const);"},
{"lineNum":"  160","line":""},
{"lineNum":"  161","line":"  KOKKOS_FUNCTION"},
{"lineNum":"  162","line":"  static task_root_type* pop_ready_task(task_root_type* volatile* const);"},
{"lineNum":"  163","line":""},
{"lineNum":"  164","line":"  KOKKOS_FUNCTION static void decrement(task_root_type* task);"},
{"lineNum":"  165","line":""},
{"lineNum":"  166","line":" public:"},
{"lineNum":"  167","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  168","line":"  int allocation_count() const noexcept { return m_count_alloc; }","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  169","line":""},
{"lineNum":"  170","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  171","line":"  void initialize_team_queues(int /*pool_size*/) const noexcept {}","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  172","line":""},
{"lineNum":"  173","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  174","line":"  task_root_type* attempt_to_steal_task() const noexcept { return nullptr; }","class":"lineNoCov","hits":"0","possible_hits":"2",},
{"lineNum":"  175","line":""},
{"lineNum":"  176","line":"  KOKKOS_INLINE_FUNCTION"},
{"lineNum":"  177","line":"  team_queue_type& get_team_queue(int /*team_rank*/) { return *this; }","class":"lineNoCov","hits":"0","possible_hits":"3",},
{"lineNum":"  178","line":""},
{"lineNum":"  179","line":"  // void execute() { specialization::execute( this ); }"},
{"lineNum":"  180","line":""},
{"lineNum":"  181","line":"  template <typename FunctorType>"},
{"lineNum":"  182","line":"  void proc_set_apply(typename task_root_type::function_type* ptr) {"},
{"lineNum":"  183","line":"    using specialization ="},
{"lineNum":"  184","line":"        TaskQueueSpecialization<BasicTaskScheduler<ExecSpace, TaskQueue>>;"},
{"lineNum":"  185","line":"    specialization::template proc_set_apply<FunctorType>(ptr);"},
{"lineNum":"  186","line":"  }"},
{"lineNum":"  187","line":""},
{"lineNum":"  188","line":"  // Assign task pointer with reference counting of assigned tasks"},
{"lineNum":"  189","line":"  KOKKOS_FUNCTION static void assign(task_root_type** const lhs,"},
{"lineNum":"  190","line":"                                     task_root_type* const rhs) {","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  191","line":"#if 0"},
{"lineNum":"  192","line":"  {"},
{"lineNum":"  193","line":"    printf( \"assign( 0x%lx { 0x%lx %d %d } , 0x%lx { 0x%lx %d %d } )\\n\""},
{"lineNum":"  194","line":"          , uintptr_t( lhs ? *lhs : 0 )"},
{"lineNum":"  195","line":"          , uintptr_t( lhs && *lhs ? (*lhs)->m_next : 0 )"},
{"lineNum":"  196","line":"          , int( lhs && *lhs ? (*lhs)->m_task_type : 0 )"},
{"lineNum":"  197","line":"          , int( lhs && *lhs ? (*lhs)->m_ref_count : 0 )"},
{"lineNum":"  198","line":"          , uintptr_t(rhs)"},
{"lineNum":"  199","line":"          , uintptr_t( rhs ? rhs->m_next : 0 )"},
{"lineNum":"  200","line":"          , int( rhs ? rhs->m_task_type : 0 )"},
{"lineNum":"  201","line":"          , int( rhs ? rhs->m_ref_count : 0 )"},
{"lineNum":"  202","line":"          );"},
{"lineNum":"  203","line":"    fflush( stdout );"},
{"lineNum":"  204","line":"  }"},
{"lineNum":"  205","line":"#endif"},
{"lineNum":"  206","line":""},
{"lineNum":"  207","line":"    if (*lhs) decrement(*lhs);","class":"lineNoCov","hits":"0","possible_hits":"3",},
{"lineNum":"  208","line":"    if (rhs) {","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  209","line":"      Kokkos::atomic_increment(&(rhs->m_ref_count));"},
{"lineNum":"  210","line":"    }"},
{"lineNum":"  211","line":""},
{"lineNum":"  212","line":"    // Force write of *lhs"},
{"lineNum":"  213","line":""},
{"lineNum":"  214","line":"    *static_cast<task_root_type* volatile*>(lhs) = rhs;","class":"lineNoCov","hits":"0","possible_hits":"4",},
{"lineNum":"  215","line":""},
{"lineNum":"  216","line":"    Kokkos::memory_fence();"},
{"lineNum":"  217","line":"  }","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  218","line":""},
{"lineNum":"  219","line":"  KOKKOS_FUNCTION"},
{"lineNum":"  220","line":"  size_t allocate_block_size(size_t n);  ///< Actual block size allocated"},
{"lineNum":"  221","line":""},
{"lineNum":"  222","line":"  KOKKOS_FUNCTION"},
{"lineNum":"  223","line":"  void* allocate(size_t n);  ///< Allocate from the memory pool"},
{"lineNum":"  224","line":""},
{"lineNum":"  225","line":"  KOKKOS_FUNCTION"},
{"lineNum":"  226","line":"  void deallocate(void* p, size_t n);  ///< Deallocate to the memory pool"},
{"lineNum":"  227","line":""},
{"lineNum":"  228","line":"  //----------------------------------------"},
{"lineNum":"  229","line":"  /**\\brief  Allocation size for a spawned task */"},
{"lineNum":"  230","line":""},
{"lineNum":"  231","line":"  template <typename FunctorType>"},
{"lineNum":"  232","line":"  KOKKOS_FUNCTION size_t spawn_allocation_size() const {"},
{"lineNum":"  233","line":"    using value_type = typename FunctorType::value_type;"},
{"lineNum":"  234","line":""},
{"lineNum":"  235","line":"    using task_type = Impl::Task<execution_space, value_type, FunctorType>;"},
{"lineNum":"  236","line":""},
{"lineNum":"  237","line":"    enum : size_t { align = (1 << 4), align_mask = align - 1 };"},
{"lineNum":"  238","line":"    enum : size_t { task_size = sizeof(task_type) };"},
{"lineNum":"  239","line":"    enum : size_t { result_size = Impl::TaskResult<value_type>::size };"},
{"lineNum":"  240","line":"    enum : size_t {"},
{"lineNum":"  241","line":"      alloc_size = ((task_size + align_mask) & ~align_mask) +"},
{"lineNum":"  242","line":"                   ((result_size + align_mask) & ~align_mask)"},
{"lineNum":"  243","line":"    };"},
{"lineNum":"  244","line":""},
{"lineNum":"  245","line":"    return m_memory.allocate_block_size(task_size);"},
{"lineNum":"  246","line":"  }"},
{"lineNum":"  247","line":""},
{"lineNum":"  248","line":"  /**\\brief  Allocation size for a when_all aggregate */"},
{"lineNum":"  249","line":""},
{"lineNum":"  250","line":"  KOKKOS_FUNCTION"},
{"lineNum":"  251","line":"  size_t when_all_allocation_size(int narg) const {","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  252","line":"    return m_memory.allocate_block_size(sizeof(task_root_type) +","class":"lineNoCov","hits":"0","possible_hits":"4",},
{"lineNum":"  253","line":"                                        narg * sizeof(task_root_type*));","class":"lineNoCov","hits":"0","possible_hits":"1",},
{"lineNum":"  254","line":"  }"},
{"lineNum":"  255","line":"};"},
{"lineNum":"  256","line":""},
{"lineNum":"  257","line":"} /* namespace Impl */"},
{"lineNum":"  258","line":"} /* namespace Kokkos */"},
{"lineNum":"  259","line":""},
{"lineNum":"  260","line":"//----------------------------------------------------------------------------"},
{"lineNum":"  261","line":"//----------------------------------------------------------------------------"},
{"lineNum":"  262","line":""},
{"lineNum":"  263","line":"#endif /* #if defined( KOKKOS_ENABLE_TASKDAG ) */"},
{"lineNum":"  264","line":"#endif /* #ifndef KOKKOS_IMPL_TASKQUEUE_HPP */"},
]};
var percent_low = 25;var percent_high = 75;
var header = { "command" : "fodo_cxx", "date" : "2021-12-15 09:03:41", "instrumented" : 13, "covered" : 0,};
var merged_data = [];
