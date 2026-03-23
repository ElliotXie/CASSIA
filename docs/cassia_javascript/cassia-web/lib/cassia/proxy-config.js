/**
 * Proxy configuration for CASSIA Web.
 *
 * Auto-detects China users (via timezone) and routes API calls through
 * a Cloudflare Worker reverse proxy. Users still BYOK.
 *
 * To disable: set localStorage key "cassia_proxy" to "off"
 * To force enable: set localStorage key "cassia_proxy" to "china"
 */

const PROXY_BASE = 'https://cassia-proxy.solodev66.workers.dev';

const PROXY_ENDPOINTS = {
  china: {
    openai: `${PROXY_BASE}/openai/v1`,
    anthropic: `${PROXY_BASE}/anthropic`,
    openrouter: `${PROXY_BASE}/openrouter/api/v1`,
  }
};

/**
 * Detect if user is likely in China based on browser timezone.
 */
function detectChina() {
  try {
    const tz = Intl.DateTimeFormat().resolvedOptions().timeZone || '';
    // Direct China timezone match
    if (tz.startsWith('Asia/Shanghai') || tz.startsWith('Asia/Chongqing') ||
        tz.startsWith('Asia/Harbin') || tz.startsWith('Asia/Urumqi') ||
        tz === 'PRC') {
      return true;
    }
    // Browser language as additional signal for UTC+8 zones
    const lang = navigator.language || '';
    if (lang.startsWith('zh') && tz.startsWith('Asia/')) {
      return true;
    }
  } catch (e) {
    // ignore
  }
  return false;
}

/**
 * Get the active proxy name, or null if disabled.
 * Checks localStorage override first, then auto-detects.
 */
export function getProxy() {
  if (typeof window !== 'undefined') {
    const override = localStorage.getItem('cassia_proxy');
    if (override === 'off' || override === 'none') return null;
    if (override && PROXY_ENDPOINTS[override]) return override;
  }
  return detectChina() ? 'china' : null;
}

/**
 * Set proxy override. Persists in localStorage.
 * @param {string|null} proxyName - "china" to enable, null/"off" to disable
 */
export function setProxy(proxyName) {
  if (typeof window !== 'undefined') {
    if (proxyName === null || proxyName === 'off') {
      localStorage.setItem('cassia_proxy', 'off');
    } else {
      localStorage.setItem('cassia_proxy', proxyName);
    }
  }
}

/**
 * Get the base URL for a provider, using proxy if active.
 * @param {string} provider - "openai", "anthropic", or "openrouter"
 * @returns {{ baseUrl: string|null, isProxied: boolean }}
 */
export function getProviderBaseUrl(provider) {
  const proxy = getProxy();
  if (proxy && PROXY_ENDPOINTS[proxy]) {
    const url = PROXY_ENDPOINTS[proxy][provider.toLowerCase()];
    if (url) return { baseUrl: url, isProxied: true };
  }
  return { baseUrl: null, isProxied: false };
}
