#include <iostream>
#include <functional>
#include <queue>
#include <string>

struct task_queue
{
    template < typename CALLABLE, typename... ARGS >
    void run_when_asked( CALLABLE fn, ARGS&&... args )
    { pending.push( std::bind( fn, args... ) ) ; }

    void execute_pending()
    {
        while( !pending.empty() )
        {
            pending.front()() ;
            pending.pop() ;
        }
    }

    std::queue< std::function< void() > > pending ;
};

int free_fun( int a, int b )
{
    std::cout << "free_fun(" << a << ',' << b << ")\n" ;
    return a+b ;
}



int main()
{
    task_queue tq ;
    tq.run_when_asked( free_fun, 1, 2 ) ;

    tq.execute_pending() ;
}
