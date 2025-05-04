document.addEventListener('DOMContentLoaded', function() {
    // Tab functionality
    const tabButtons = document.querySelectorAll('.tab-btn');
    const tabContents = document.querySelectorAll('.tab-content');

    tabButtons.forEach(button => {
        button.addEventListener('click', () => {
            // Remove active class from all buttons and contents
            tabButtons.forEach(btn => btn.classList.remove('active'));
            tabContents.forEach(content => content.classList.remove('active'));
            
            // Add active class to clicked button and corresponding content
            button.classList.add('active');
            const tabId = button.getAttribute('data-tab');
            document.getElementById(tabId).classList.add('active');
        });
    });

    // Language switching
    const enSwitch = document.getElementById('en-switch');
    const zhSwitch = document.getElementById('zh-switch');

    if (enSwitch && zhSwitch) {
        enSwitch.addEventListener('click', (e) => {
            e.preventDefault();
            if (!enSwitch.classList.contains('active')) {
                enSwitch.classList.add('active');
                zhSwitch.classList.remove('active');
                // Here you would implement the actual language switch
                // For now, we'll just show an alert
                window.location.href = 'index.html';
            }
        });

        zhSwitch.addEventListener('click', (e) => {
            e.preventDefault();
            if (!zhSwitch.classList.contains('active')) {
                zhSwitch.classList.add('active');
                enSwitch.classList.remove('active');
                // Here you would implement the actual language switch
                // For now, we'll just show an alert
                window.location.href = 'index_cn.html';
            }
        });
    }

    // Smooth scrolling for anchor links
    document.querySelectorAll('a[href^="#"]').forEach(anchor => {
        anchor.addEventListener('click', function(e) {
            if (this.getAttribute('href') === '#') return;
            
            e.preventDefault();
            const targetId = this.getAttribute('href');
            
            if (targetId !== '#') {
                const targetElement = document.querySelector(targetId);
                if (targetElement) {
                    window.scrollTo({
                        top: targetElement.offsetTop - 80, // Offset for header
                        behavior: 'smooth'
                    });
                }
            }
        });
    });

    // Add scroll animation for elements
    const fadeInElements = document.querySelectorAll('.card, .news-item, .model-card');
    
    const fadeInOptions = {
        threshold: 0.1,
        rootMargin: "0px 0px -100px 0px"
    };
    
    const fadeInObserver = new IntersectionObserver(function(entries, observer) {
        entries.forEach(entry => {
            if (entry.isIntersecting) {
                entry.target.style.opacity = 1;
                entry.target.style.transform = 'translateY(0)';
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
}); 