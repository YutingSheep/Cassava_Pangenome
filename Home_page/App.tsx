import { Header } from './components/Header';
import { Hero } from './components/Hero';
import { ModulesSection } from './components/ModulesSection';
import { Footer } from './components/Footer';
import { useI18nStore } from './lib/i18n';
import './App.css';

function App() {
  const { language } = useI18nStore();

  return (
    <div className={`App ${language}`}>
      <Header />
      <main>
        <Hero />
        <ModulesSection />
      </main>
      <Footer />
    </div>
  );
}

export default App;
