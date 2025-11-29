document.addEventListener('DOMContentLoaded', function() {
    // Tab functionality
    const tabsContainer = document.querySelector('.tabs');
    if (tabsContainer) {
        tabsContainer.addEventListener('click', function(e) {
            if (e.target.classList.contains('tab-btn')) {
                const button = e.target;
                const tabButtons = document.querySelectorAll('.tab-btn');
                const tabContents = document.querySelectorAll('.tab-content');
                
                // Remove active class from all buttons and contents
                tabButtons.forEach(btn => btn.classList.remove('active'));
                tabContents.forEach(content => content.classList.remove('active'));
                
                // Add active class to clicked button and corresponding content
                button.classList.add('active');
                const tabId = button.getAttribute('data-tab');
                document.getElementById(tabId).classList.add('active');
            }
        });
    }

    // Language switching
    const languageSwitcher = document.querySelector('.language-switch');
    if (languageSwitcher) {
        languageSwitcher.addEventListener('click', function(e) {
            e.preventDefault();
            
            if (e.target.id === 'en-switch' && !e.target.classList.contains('active')) {
                window.location.href = 'index.html';
            } else if (e.target.id === 'zh-switch' && !e.target.classList.contains('active')) {
                window.location.href = 'index_cn.html';
            }
        });
    }

    // Smooth scrolling for anchor links with performance optimization
    document.addEventListener('click', function(e) {
        const target = e.target.closest('a[href^="#"]');
        if (!target) return;
        
        if (target.getAttribute('href') === '#') return;
        
        e.preventDefault();
        const targetId = target.getAttribute('href');
        
        if (targetId !== '#') {
            const targetElement = document.querySelector(targetId);
            if (targetElement) {
                window.scrollTo({
                    top: targetElement.offsetTop - 80, // Offset for header
                    behavior: 'smooth'
                });
            }
        }
    }, { passive: false });

    // Add scroll animation for elements using Intersection Observer API
    const fadeInElements = document.querySelectorAll('.card, .news-item, .model-card');
    
    if ('IntersectionObserver' in window && fadeInElements.length > 0) {
        const fadeInOptions = {
            threshold: 0.1,
            rootMargin: "0px 0px -100px 0px"
        };
        
        const fadeInObserver = new IntersectionObserver(function(entries, observer) {
            entries.forEach(entry => {
                if (entry.isIntersecting) {
                    // Use requestAnimationFrame for smoother animations
                    requestAnimationFrame(() => {
                        entry.target.style.opacity = 1;
                        entry.target.style.transform = 'translateY(0)';
                    });
                    observer.unobserve(entry.target);
                }
            });
        }, fadeInOptions);
        
        fadeInElements.forEach(element => {
            element.style.opacity = 0;
            element.style.transform = 'translateY(20px)';
            element.style.transition = 'opacity 0.5s ease, transform 0.5s ease';
            fadeInObserver.observe(element);
        });
    }
}); 